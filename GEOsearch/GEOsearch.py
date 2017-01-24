import csv, os, re, psutil, urllib, contextlib, gc
from collections import defaultdict
from GSM import GSM, search_term_to_GSM
from input_search_utils import SOFTQuickRelated, input_finder, has_features
from pickleUtils import load_obj
from multiprocessing import Process, Queue


def get_MetaInfo(path):
    proc = psutil.Process()
    file_obj = open(path, "r")
    info = file_obj.readlines()
    file_obj.close()
    if len(proc.open_files()) > 3:
        return None
    del file_obj
    # if count % 50 == 0:
    #     gc.collect()
    return info


def get_WebInfo(url):
    with contextlib.closing(urllib.urlopen(url)) as web:
        info = web.readlines()
    web.close()
    del web
    # if count % 50 == 0:
    #     gc.collect()
    return info


def SOFTQuickParser(output_surfix, features, features_begin,
                    type_seq="chip-seq", ignorecase=True, geo=False, geofile=None, output_type="Homo sapiens",
                    encode_remove=True, roadmap_remove=True, encode_pkl=None, roadmap_pkl=None, GSMGSE_pkl=None,
                    cwd=None, process=20):

    encodeGSE = load_obj(encode_pkl)

    roadmapGSE = load_obj(roadmap_pkl)

    excludedGSE = set()
    if encode_remove:
        excludedGSE = excludedGSE.union(encodeGSE)
    if roadmap_remove:
        excludedGSE = excludedGSE.union(roadmapGSE)

    GSEGSM_map = load_obj(GSMGSE_pkl)

    excludedGSM = set()

    for gse in excludedGSE:
        excludedGSM = excludedGSM.union(GSEGSM_map[gse])

    type_seq = type_seq.lower()
    # print len(map)

    geoGSMs = set()
    if geofile is not None:
        file_obj = open(geofile, "r")
        for line in file_obj.readlines():
            geoGSMs.add(line.strip())
        file_obj.close()
    else:
        geoGSMs = search_term_to_GSM(features+features_begin)

        if (cwd is not None) and (not geo):
            localGSMs = set([x[:-4] for x in os.listdir(cwd) if x.startswith("GSM")])
            geoGSMs = geoGSMs.union(localGSMs)
    geoGSMs = list(geoGSMs)

    print "total ", len(geoGSMs), " candidate found in search", features

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
    notFeature = {}
    for i in range(process):
        cur_samples, cur_Human_Samples, cur_notFeature = queue.get()
        samples.update(cur_samples)
        Human_Samples.update(cur_Human_Samples)
        notFeature.update(cur_notFeature)
    for p in processes:
        p.join()

    print "total ", output_type, " sample found", len(Human_Samples)

    groupByGSE, encodeGSE, relatedSamples = SOFTQuickRelated(Human_Samples, output_type, type_seq,
                                                             GSEGSM_map, encode_remove, encodeGSE, cwd, process)

    first_category, third_category = input_finder(output_surfix, Human_Samples, groupByGSE, encodeGSE, relatedSamples,
                                                  features, features_begin, ignorecase, output_type)

    #### output results to csv
    output_type = output_type.replace(" ", "_")

    outputHuman = "./"+"GEOsearch"+ output_type +"With" + output_surfix + ".csv"
    outputSample = "./"+"GEOsearch"+"sampleWith" + output_surfix + ".csv"
    outputNoFeature = "./"+"GEOsearch"+"noWith" + output_surfix+ ".csv"

    csv_file = open(outputSample, "wb")
    writer = csv.writer(csv_file)
    writer.writerow(
        ['Sample_ID', "Series_ID", output_surfix + "_description", "Organism", "Title", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
         "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab",
         "Feature in Title or Ab"])
    for sample in samples.values():
        writer.writerow(
            [sample.id, sample.series, sample.features, sample.organism, sample.title, sample.platForm, sample.InstrumentID,
             sample.SRA, sample.libraryStrategy, sample.tissue, sample.cellLine, sample.cellType,
             sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found,
             sample.title_ab])
    csv_file.close()


    csv_file = open(outputNoFeature, "wb")
    writer = csv.writer(csv_file)
    writer.writerow(
        ['Sample_ID', "Series_ID", output_surfix + "_description", "Organism", "Title", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
         "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab",
         "Feature in Title or Ab"])
    for sample in notFeature.values():
        writer.writerow(
            [sample.id, sample.series, sample.features, sample.organism, sample.title, sample.platForm, sample.InstrumentID,
             sample.SRA, sample.libraryStrategy, sample.tissue, sample.cellLine, sample.cellType,
             sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found,
             sample.title_ab])
    csv_file.close()


    csv_file = open(outputHuman, "wb")
    writer = csv.writer(csv_file)
    writer.writerow(
        ['Sample_ID', "Series_ID", output_surfix + "_description", "Input_ID", "Input_Description", "Organism", "Title",
         "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
         "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title",
         "Feature in Ab", "Feature in Title or Ab"])

    for sample in Human_Samples.values():
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

        writer.writerow(
            [sample.id, sample.series, sample.features, potential_input_id, potential_input_title, sample.organism, sample.title,
             sample.platForm, sample.InstrumentID, sample.SRA, sample.libraryStrategy, sample.tissue,
             sample.cellLine, sample.cellType, sample.disease, sample.treatment, sample.genotype, sample.antibody,
             sample.title_found, sample.ab_found, sample.title_ab])
    csv_file.close()

    return Human_Samples


def feature_filter(geoGSMs, queue, features, features_begin, excludedGSM,
                    type_seq, ignorecase, output_type,
                    cwd):
    print os.getpid()
    samples = {}
    Human_Samples = {}
    notFeature = {}

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

        if cwd is not None:
            if not cwd.endswith("/"):
                cwd += "/"

            if os.path.isfile(cwd + sampleName + ".txt"):
                info = get_MetaInfo(cwd + sampleName + ".txt")
            else:
                info = get_WebInfo(sample.url)
        else:
            info = get_WebInfo(sample.url)

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

            if key.lower() in ["tissue", "body part", "body site"]:
                tissue = value

            if key.lower() in ["cancer type", "tumor type", "tumor region", "disease", "disease state",
                               "disease status"]:
                disease = value

            if key.lower() in ["background strain", "strain", "strain number", "mouse strain", "strain background",
                               "cell line background", "genetic background", "genotype", "genotype/variation",
                               "strain/background", "variation"]:
                genoType[key] = value

            if key.lower() in ["cell line", "cell", "cells pointed by barcodes",
                               "chicken line", "line"]:
                cellLine += value

            if key.lower() in ["cell_type", "cell-type", "cell type", "cell lineage"]:
                cellType += value

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
                and sample.id not in excludedGSM:
            if sample.title_ab:
                if sample.title.lower().find("input") == -1 \
                        and sample.title.lower().find("wce") == -1 \
                        and sample.title.find("IgG") == -1:
                    Human_Samples[sample.id] = sample
                else:
                    print sample.title, "title has input or wce or IgG!"
                    continue

        if len(target_feature) != 0 and sample.id not in excludedGSM:
            samples[sampleName] = sample
        else:
            notFeature[sampleName] = sample
    queue.put((samples, Human_Samples, notFeature))



if __name__ == "__main__":
    SOFTQuickParser("H3K27me3", ["H3K27me3"], [], geo=True, geofile='./testlist.txt', encode_remove=True,encode_pkl='./pkl/ENCODE_gse.pkl', roadmap_pkl='./pkl/Roadmap_gse.pkl', GSMGSE_pkl='./pkl/GSMGSE_map.pkl')