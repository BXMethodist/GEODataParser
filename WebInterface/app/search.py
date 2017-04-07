"""
Copyright (c) <2017> <Dr. Kaifu Chen lab, Research Institute, Houston Methodist Hospital >

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import re, sqlite3, pandas as pd
from collections import defaultdict
from sample import GSM, search_term_to_GSM
from input_search_utils import SOFTQuickRelated, input_finder, has_features, get_MetaInfo, get_WebInfo
from pickleUtils import load_obj
from multiprocessing import Process, Queue
from encode import encode_search
import re


def SOFTQuickParser(output_surfix, output_path, features, features_begin,
                    type_seq="chip-seq", ignorecase=True, geo=False, geofile=None, output_type="Homo sapiens",
                    encode_remove=True, roadmap_remove=True, encode_pkl=None, roadmap_pkl=None, GGRmap_pkl=None,
                    GSMGSE_pkl=None, cwd=None, process=20, email=None):

    encodeGSE = load_obj(encode_pkl)

    roadmapGSE = load_obj(roadmap_pkl)

    GGRmapGSE = load_obj(GGRmap_pkl)

    excludedGSE = set()
    if encode_remove:
        excludedGSE = excludedGSE.union(encodeGSE)
        excludedGSE = excludedGSE.union(roadmapGSE)
        excludedGSE = excludedGSE.union(GGRmapGSE)

    GSEGSM_map = load_obj(GSMGSE_pkl)

    excludedGSM = set()

    for gse in excludedGSE:
        excludedGSM = excludedGSM.union(GSEGSM_map[gse])

    type_seq = type_seq.lower()
    # print len(map)

    geoGSMs = set()
    encode_candidates = set()
    if geofile is not None:
        file_obj = open(geofile, "r")
        for line in file_obj.readlines():
            if line.startswith('ENC'):
                encode_candidates.add(line.strip())
            else:
                geoGSMs.add(line.strip())
        file_obj.close()
    else:
        geoGSMs = search_term_to_GSM(features+features_begin)

        if (cwd is not None) and (not geo):
            db = sqlite3.connect(cwd)
            db.text_factory = str
            query = db.execute('SELECT GSM_ID from GSM').fetchall()
            localGSMs = set([x[0] for x in query])
            geoGSMs = geoGSMs.union(localGSMs)
    geoGSMs = list(geoGSMs)

    # print "total ", len(geoGSMs), " candidate found in search", features
    ## test only
    # geoGSMs = []

    ####

    queue = Queue()
    chucksize = len(geoGSMs)/(process-1) if process > 1 else len(geoGSMs)
    processes = []
    for i in range(process):
        cur_geoGSMs = geoGSMs[i*chucksize:(i+1)*chucksize]
        p = Process(target=feature_filter, args=(cur_geoGSMs, queue, features, features_begin, excludedGSM,
                    type_seq, ignorecase, output_type,cwd,))
        processes.append(p)
        p.start()

    samples = {}
    Human_Samples = {}

    for i in range(process):
        cur_samples, cur_Human_Samples = queue.get()
        samples.update(cur_samples)
        Human_Samples.update(cur_Human_Samples)
    for p in processes:
        p.join()
    # print "total ", output_type, " sample found", len(Human_Samples)

    # looking for input
    if type_seq == 'chip-seq':
        groupByGSE, excludedGSE, relatedSamples = SOFTQuickRelated(Human_Samples, output_type, type_seq,
                                                                 GSEGSM_map, encode_remove, excludedGSE, cwd, process)

        first_category, third_category = input_finder(output_surfix, output_path, Human_Samples, groupByGSE, excludedGSE, relatedSamples,
                                                      features, features_begin, ignorecase, output_type)
    else:
        first_category, third_category = defaultdict(set), defaultdict(set)
    #

    # search data in encode database
    if geofile is not None and len(encode_candidates) != 0:
        samples_encode, human_encode, human_encode_map = encode_search(output_surfix, features,
                                                                       keywords_begin=features_begin,
                                                                       type_seq=type_seq,
                                                                       candidates=encode_candidates,
                                                                       ignorecase=ignorecase,
                                                                       output_type=output_type)
    elif geofile is not None and len(encode_candidates) == 0:
        samples_encode, human_encode, human_encode_map = encode_search(output_surfix, features,
                                                                       keywords_begin=features_begin,
                                                                       type_seq=type_seq,
                                                                       candidates=None,
                                                                       ignorecase=ignorecase,
                                                                       output_type=output_type)
    elif geofile is None:
        samples_encode, human_encode, human_encode_map = encode_search(output_surfix, features,
                                                                       keywords_begin=features_begin,
                                                                       type_seq=type_seq,
                                                                       candidates=encode_candidates,
                                                                       ignorecase=ignorecase,
                                                                       output_type=output_type)

    # output results to csv
    output_type = output_type.replace(" ", "_")

    if not output_path.endswith("/"):
        output_path += "/"
    outputHuman = output_path+"Search_Result" + output_type + "With" + output_surfix + ".csv"
    outputSample = output_path+"Search_Result" + "sampleWith" + output_surfix + ".csv"

    table = []
    headers = ['Data_ID', "Study_ID", "Data_Description", "Title",
               "Instrument_Model", "Raw Data", "Sequencing_Protocol", "Species", "Cell Line", "Cell Type",
               "Experiment target/antibody", 'Confidence', 'Tissue', 'Organ']

    for sample in samples.values():
        if sample.title_found and sample.ab_found:
            confidence = 'High Confident'
        elif sample.title_found:
            confidence = 'Medium Confident'
        elif sample.ab_found:
            confidence = 'Low Confident'
        else:
            confidence = 'No Confident'

        features = str(sample.features)[1:-1].replace("u'", "")
        features = features.replace("'", '')
        antibody = str(sample.antibody)[1:-1].replace("u'", "")
        antibody = antibody.replace("'", '')

        row = [sample.id, ",".join(list(sample.series)),
               features,
               sample.title,
               sample.InstrumentID, sample.SRA, sample.libraryStrategy, sample.organism, sample.cellLine, sample.cellType,
               antibody, confidence, sample.tissue, sample.organ]
        table.append(row)

    df = pd.DataFrame(table, columns=headers)
    df = df.set_index(['Data_ID'])
    if samples_encode is not None:
        df = df.append(samples_encode)
    df.to_csv(outputSample, sep=',', encoding='utf-8')

    table = []
    headers = ['Data_ID', "Study_ID", "Data_Description", "Title",
               "Input", "Input_Description", "Instrument_Model", "Raw Data", "Sequencing_Protocol",
               "Species", "Cell Line", "Cell Type",
               "Experiment target/antibody", 'Confidence', 'Tissue', 'Organ']

    for sample in Human_Samples.values():
        if sample.title_found and sample.ab_found:
            confidence = 'High Confident'
        elif sample.title_found:
            confidence = 'Medium Confident'
        elif sample.ab_found:
            confidence = 'Low Confident'
        else:
            confidence = 'No Confident'

        features = str(sample.features)[1:-1].replace("u'", "")
        features = features.replace("'", '')
        antibody = str(sample.antibody)[1:-1].replace("u'", "")
        antibody = antibody.replace("'", '')

        potential_input_id = ""
        potential_input_title = ""
        if len(first_category[sample.id]) != 0:
            for id in first_category[sample.id]:
                potential_input_id += id + ","
                potential_input_title += relatedSamples[id].title +","
            potential_input_id = potential_input_id[:-1]
            potential_input_title = potential_input_title[:-1]

        elif len(third_category) != 0:
            for id in third_category[sample.id]:
                potential_input_id += id + ","
                for key, value in relatedSamples[id].antibody.items():
                    potential_input_title += key+":"+value+","
            potential_input_id = potential_input_id[:-1]
            potential_input_title = potential_input_title[:-1]

        row = [sample.id, ",".join(list(sample.series)),
               features,
               sample.title,
               potential_input_id, potential_input_title, sample.InstrumentID, sample.SRA, sample.libraryStrategy,
               sample.organism, sample.cellLine, sample.cellType,
               antibody, confidence, sample.tissue, sample.organ]
        table.append(row)

    df = pd.DataFrame(table, columns=headers)
    df = df.set_index(['Data_ID'])

    if human_encode is not None:
        df = df.append(human_encode)

    high_human = output_path + "Search_Result" + output_type + "With" + output_surfix + '_High_Confidence'+ ".csv"
    high_human_df = df[df['Confidence'] == 'High Confident']
    medium_human = output_path + "Search_Result" + output_type + "With" + output_surfix + '_Medium_Confidence' + ".csv"
    medium_human_df = df[df['Confidence'] == 'Medium Confident']
    low_human = output_path + "Search_Result" + output_type + "With" + output_surfix + '_Low_Confidence' + ".csv"
    low_human_df = df[df['Confidence'] == 'Low Confident']
    no_human = output_path + "Search_Result" + output_type + "With" + output_surfix + '_No_Confidence' + ".csv"
    no_human_df = df[df['Confidence'] == 'No Confident']

    df.to_csv(outputHuman, sep=',', encoding='utf-8')
    high_human_df.to_csv(high_human, sep=',', encoding='utf-8')
    medium_human_df.to_csv(medium_human, sep=',', encoding='utf-8')
    low_human_df.to_csv(low_human, sep=',', encoding='utf-8')
    no_human_df.to_csv(no_human, sep=',', encoding='utf-8')

    if human_encode_map is not None:
        Human_Samples.update(human_encode_map)
    return Human_Samples


def feature_filter(geoGSMs, queue, features, features_begin, excludedGSM,
                    type_seq, ignorecase, output_type, cwd):
    # print "Process id is ", os.getpid()
    samples = {}
    Human_Samples = {}

    if cwd is None:
        db = None
    else:
        db = sqlite3.connect(cwd)
        # db.text_factory = str

    count = 0

    for sampleName in geoGSMs:
        sample = GSM(sampleName)

        sampleTitle = ""
        sampleType = ""
        sampleLibraryStrategy = ''
        sampleOrganism = ''
        samplePlatForm = ''
        sampleInstrumentID = ''

        characteristics = defaultdict(str)
        supplementaryData = defaultdict(str)
        relations = defaultdict(str)
        sampleSeriesID = set()
        target_feature = {}

        antibody = {}
        treatment = {}
        tissue = ""
        disease = ""
        cellLine = ""
        cellType = ""
        genoType = {}
        title_found = False
        ab_found = False

        info = get_MetaInfo(db, sample, count)

        count += 1

        if info is None:
            continue

        for line in info:
            if line.startswith("!Sample_title"):
                sampleTitle = line[line.find("=") + 1:].strip()
                if sampleTitle.find(";") != -1:
                    sampleTitle = sampleTitle[:sampleTitle.find(";")]
                if has_features(sampleTitle, features, features_begin, ignorecase):
                    target_feature["Title"] = sampleTitle
                    title_found = True
            if line.startswith("!Sample_type"):
                sampleType = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_organism"):
                sampleOrganism = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_characteristics_ch"):
                characteristic = line[line.find("=") + 1:].strip()
                key, value = characteristic[:characteristic.find(":")].strip(), characteristic[
                                                                                characteristic.find(":") + 1:].strip()
                if key in characteristics:
                    characteristics[key] += ", " + value
                else:
                    characteristics[key] = value
                if has_features(value, features, features_begin, ignorecase):
                    target_feature[key] = value
            if line.startswith("!Sample_platform_id "):
                samplePlatForm = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_library_strategy"):
                sampleLibraryStrategy = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_supplementary"):
                dataUrl = line[line.find("=") + 1:].strip()
                keyName = "Data" + str(len(supplementaryData) + 1)
                supplementaryData[keyName] = dataUrl
            if line.startswith("!Sample_relation"):
                relation = line[line.find("=") + 1:].strip()
                key, value = relation[:relation.find(":")].strip(), relation[relation.find(":") + 1:].strip()
                relations[key] = value
            if line.startswith("!Sample_series_id"):
                sampleSeriesID.add(line[line.find("=") + 1:].strip())
            if line.startswith("!Sample_instrument_model"):
                sampleInstrumentID = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_contact_institute = ENCODE DCC"):
                sample.encode = True

        sample.characteristics = characteristics
        sample.supplementaryData = supplementaryData
        sample.title = sampleTitle
        sample.type = sampleType
        sample.libraryStrategy = sampleLibraryStrategy
        sample.organism = sampleOrganism
        sample.SRA = relations["SRA"]
        sample.series = list(sampleSeriesID)

        sample.platForm = samplePlatForm
        sample.features = target_feature
        sample.InstrumentID = sampleInstrumentID

        for key, value in characteristics.items():
            if key.lower() in ["chip antibody", "chip", "antigen", "antibody", "antibodies", "chip antibodies",
                               "antibody name", "antibody target", "target", "antibody/capture",
                               "antibody/vendor/catalog#",
                               "chip ab", "chip antibody1", "chip antibody2", "chip-antibody", "chip_antibody",
                               "chip-seq antibody", "chip_antibodies", "chip-antibodies", "histone mark",
                               "epigenetic feature",
                               "histone modification", "antibody antibodydescription", "chip antibody (epitope/name)",
                               "factor", "chip antibody/mbd affinity column", "chip/dip antibody", "antibody epiptope",
                               "antibody source", 'modification', "antibody (vendor': ' catalog#, or reference)",
                               "experiment", "purification antibody", "antibody/details", "antibody epiptope",
                               "antibody information", "chip antibody / digestive enzyme", "chip antiboy",
                               "ip antibody", "chip antibody target", "modification", "histone", "enrichment procedure",
                               "antibody (vendor': ' catalog#, or reference)",
                               "developmental stage/condition/extract protocol",
                               "antibody source"] \
                    or re.search('antibody epitope|using[\w\s]+antibod|immunoprecipitat', key, flags=re.IGNORECASE):
                if key in antibody:
                    antibody[key] += ", " + value
                else:
                    antibody[key] = value
            if key.lower() in ["treatment", "condition", "activation stimuli", "cell condition", "cell treatment",
                               "cell-treatment", "drug treatment", "stress", "overexpression", "treatment drug",
                               "treatment group"] \
                    or re.search("(?:dsrna|infect|rnai|shrna|sirna|transduc|transfec|agent[#]*[0-9]*|activat)", key,
                                 flags=re.IGNORECASE):
                treatment[key] = value

            if key.lower() in ["tissue", "body part", "body site", "tissue type"]:
                tissue = value

            if key.lower() in ["cancer type", "tumor type", "tumor region", "disease", "disease state",
                               "disease status"]:
                disease = value

            if key.lower() in ["background strain", "strain", "strain number", "mouse strain", "strain background",
                               "genetic background", "genotype", "genotype/variation",
                               "strain/background", "variation"]:
                genoType[key] = value

            if key.lower() in ["cell line", "cell", "cells pointed by barcodes", "cell line background",
                               "chicken line", "line"]:
                cellLine += value

            if key.lower() in ["cell_type", "cell-type", "cell type", "cell lineage"]:
                cellType += value

        organ_list = ['Bone', 'Cartilage', 'Fibrous joint', 'Cartilaginous joint', 'Synovial joint', 'Muscle', 'Tendon', 'Diaphragm', 'Heart', 'Bone marrow', 'Thymus', 'Spleen', 'Lymph node', 'CNS', 'Brain', 'Spinal cord', 'Ear', 'Eye', 'Integumentary system', 'Skin', 'Subcutaneous tissue', 'Breast', 'Mammary gland', 'Myeloid', 'Lymphoid', 'Nose', 'Nasopharynx', 'Larynx', 'Trachea', 'Bronchus', 'Lung', 'Mouth', 'Salivary', 'Tongue', 'Oropharynx', 'Laryngopharynx', 'Esophagus', 'Stomach', 'intestine', 'Appendix', 'Colon', 'Rectum', 'Anus', 'Liver', 'Biliary tract', 'Pancreas', 'Genitourinary', 'Kidney', 'Ureter', 'Bladder', 'Urethra', 'Uterus', 'Vagina', 'Vulva', 'Ovary', 'Placenta', 'Scrotum', 'Penis', 'Prostate', 'Testicle', 'Seminal vesicle', 'Pituitary', 'Pineal', 'Thyroid', 'Parathyroid', 'Adrenal', 'Islets of Langerhans']

        for organ in organ_list:
            if cellLine.lower().find(organ.lower()) != -1:
                sample.organ = organ
                break
            if cellType.lower().find(organ.lower()) != -1:
                sample.organ = organ
                break
            if tissue.lower().find(organ.lower()) != -1:
                sample.organ = organ
                break

        sample.antibody = antibody
        for value in sample.antibody.values():
            if has_features(value, features, features_begin, ignorecase):
                ab_found = True
                break

        sample.treatment = treatment
        sample.disease = disease
        sample.cellLine = cellLine
        sample.cellType = cellType
        sample.genotype = genoType
        sample.tissue = tissue
        sample.title_found = title_found
        sample.ab_found = ab_found
        if title_found or ab_found:
            sample.title_ab = True

        if (sample.organism.lower() == output_type.lower() or output_type is None) and (
                sample.SRA != None and sample.SRA.strip() != "") and \
                sample.InstrumentID.startswith('Illu') and (
                sample.libraryStrategy.lower() == type_seq or type_seq is None) \
                and sample.id not in excludedGSM and not sample.encode:
            if sample.title_ab:
                if sample.title.lower().find("input") == -1 \
                        and sample.title.lower().find("wce") == -1 \
                        and sample.title.find("IgG") == -1:
                    Human_Samples[sample.id] = sample
                elif type_seq.lower() != 'chip-seq':
                    Human_Samples[sample.id] = sample
                else:
                    print sample.title, "title has input or wce or IgG!"
                    continue

        if sample.title_ab:
            if sample.title.lower().find("input") != -1 \
                    or sample.title.lower().find("wce") != -1 \
                    or sample.title.find("IgG") != -1:

                sample.title_found = False
                sample.ab_found = False
                sample.title_ab = False

        if sample.title_ab and sample.id not in excludedGSM and not sample.encode:
            samples[sampleName] = sample
        # elif sample.id not in excludedGSM:
        #     print sample.id, sample.title

    # samples = {}
    queue.put((samples, Human_Samples))
    if db is not None:
        db.close()
    return
