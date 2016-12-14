from GSM import GSM
import os
import re
from collections import defaultdict
import csv
from difflib import SequenceMatcher
import sqlite3

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def SOFTQuickRelated(cwd=None, geo=True):
    if cwd == None:
        return

    Human_Samples = {}

    relatedSamples = {}

    groupByGSE = defaultdict(set)

    featureGSMs = set()

    if geo:
        file = open("./GEOsearchhumanWithH3K4me3.csv", "r")
    else:
        file = open("./humanWithH3K4me3.csv", "r")

    for line in file.readlines():
        gsmid = line.split(",")[0]
        if gsmid.startswith("GSM"):
            featureGSMs.add(gsmid)
    file.close()

    db = sqlite3.connect('/home/tmhbxx3/archive/GEO_MetaDatabase/geoMetaData.db')
    db.text_factory = str

    relatedGSEs = set()

    for gsmid in featureGSMs:
        query = db.execute('select GSE_ID from GSEtoGSM where GSM_ID = "'+gsmid+'"').fetchall()
        for gseid in query:
            relatedGSEs.add(gseid[0])

    encodeGSE = set()
    allGSEs = db.execute("select distinct GSE_ID, organization from GSE").fetchall()
    for gse in allGSEs:
        gseid, organization = gse
        if organization.find("ENCODE") != -1:
            encodeGSE.add(gseid)

    allrelatedGSEs = encodeGSE.union(relatedGSEs)

    allrelatedGSMs = set()

    for gseid in allrelatedGSEs:
        query = db.execute('select GSM_ID from GSEtoGSM where GSE_ID = "'+gseid+'"').fetchall()
        for gsmid in query:
            allrelatedGSMs.add(gsmid[0])
    db.close()

    print len(allrelatedGSMs)

    for filename in os.listdir(cwd):
        if not filename.startswith("GSM") or filename not in allrelatedGSMs:
            continue


        file = open(cwd+'/'+filename, "r")
        characteristics = defaultdict(str)
        supplementaryData = defaultdict(str)
        relations = defaultdict(str)
        sampleSeriesID = set()
        feature = {}

        antibody = {}
        treatment = {}
        tissue = None
        disease = None
        cellLine = ""
        cellType = ""
        genoType = {}
        title_found = False
        ab_found = False

        for line in file.readlines():
            if line.startswith("^SAMPLE"):
                sampleName = line[line.find("=")+1:].strip()
            if line.startswith("!Sample_title"):
                sampleTitle = line[line.find("=")+1:].strip()
                if sampleTitle.find(";"):
                    sampleTitle = sampleTitle[:sampleTitle.find(";")]
                if re.search("h3k4me3", sampleTitle, flags=re.IGNORECASE) or re.search("k4me3", sampleTitle, flags=re.IGNORECASE):
                    feature["Title"] = sampleTitle
                    title_found = True
            if line.startswith("!Sample_type"):
                sampleType = line[line.find("=")+1:].strip()
            if line.startswith("!Sample_organism"):
                sampleOrganism = line[line.find("=")+1:].strip()
            if line.startswith("!Sample_characteristics_ch"):
                characteristic = line[line.find("=")+1:].strip()
                key, value = characteristic[:characteristic.find(":")].strip(), characteristic[characteristic.find(":")+1:].strip()
                if key in characteristics:
                    characteristics[key] += ", " + value
                else:
                    characteristics[key] = value
                if re.search("h3k4me3", value, flags=re.IGNORECASE) or re.search("k4me3", value, flags=re.IGNORECASE):
                    feature[key] = value
            if line.startswith("!Sample_platform_id "):
                samplePlatForm = line[line.find("=")+1:].strip()
            if line.startswith("!Sample_library_strategy"):
                sampleLibraryStrategy = line[line.find("=")+1:].strip()
            if line.startswith("!Sample_supplementary"):
                dataUrl = line[line.find("=")+1:].strip()
                keyName = "Data" + str(len(supplementaryData)+1)
                supplementaryData[keyName] = dataUrl
            if line.startswith("!Sample_relation"):
                relation = line[line.find("=")+1:].strip()
                key, value = relation[:relation.find(":")].strip(), relation[relation.find(":") + 1:].strip()
                relations[key] = value
            if line.startswith("!Sample_series_id"):
                sampleSeriesID.add(line[line.find("=")+1:].strip())
            if line.startswith("!Sample_instrument_model"):
                sampleInstrumentID = line[line.find("=")+1:].strip()

        sample = GSM(sampleName)
        sample.characteristics = characteristics
        sample.supplementaryData = supplementaryData
        sample.title = sampleTitle
        sample.type = sampleType
        sample.libraryStrategy = sampleLibraryStrategy
        sample.organism = sampleOrganism
        sample.SRA = relations["SRA"]
        sample.series = list(sampleSeriesID)

        sample.platForm = samplePlatForm
        sample.features = feature
        sample.InstrumentID = sampleInstrumentID

        for key, value in characteristics.items():
            if key.lower() in ["chip antibody", "chip", "antigen", "antibody", "antibodies", "chip antibodies",
                               "antibody name", "antibody target", "target", "antibody/capture", "antibody/vendor/catalog#",
                               "chip ab", "chip antibody1", "chip antibody2", "chip-antibody", "chip_antibody",
                               "chip-seq antibody", "chip_antibodies", "chip-antibodies", "histone mark", "epigenetic feature",
                               "histone modification", "antibody antibodydescription", "chip antibody (epitope/name)",
                               "factor", "chip antibody/mbd affinity column", "chip/dip antibody", "antibody epiptope",
                               "antibody source", 'modification', "antibody (vendor': ' catalog#, or reference)",
                               "experiment", "purification antibody", "antibody/details", "antibody epiptope",
                               "antibody information", "chip antibody / digestive enzyme", "chip antiboy",
                               "ip antibody", "chip antibody target", "modification", "histone", "enrichment procedure",
                               "antibody (vendor': ' catalog#, or reference)", "developmental stage/condition/extract protocol",
                               "antibody source"] \
                    or re.search('antibody epitope|using[\w\s]+antibod|immunoprecipitat', key, flags=re.IGNORECASE):
                if key in antibody:
                    antibody[key] += ", " + value
                else:
                    antibody[key] = value
            if key.lower() in ["treatment", "condition", "activation stimuli", "cell condition", "cell treatment",
                               "cell-treatment", "drug treatment", "stress", "overexpression", "treatment drug",
                               "treatment group"] \
                    or re.search("(?:dsrna|infect|rnai|shrna|sirna|transduc|transfec|agent[#]*[0-9]*|activat)", key, flags=re.IGNORECASE):
                treatment[key] = value

            if key.lower() in ["tissue", "body part", "body site"]:
                tissue = value

            if key.lower() in ["cancer type", "tumor type", "tumor region", "disease", "disease state", "disease status"]:
                disease = value

            if key.lower() in ["background strain", "strain", "strain number", "mouse strain", "strain background",
                               "cell line background", "genetic background", "genotype", "genotype/variation",
                               "strain/background", "variation"]:
                genoType[key] = value

            if key.lower() in ["cell line",  "cell", "cells pointed by barcodes",
                           "chicken line", "line"]:
                cellLine += value

            if key.lower() in ["cell_type", "cell-type", "cell type", "cell lineage"]:
                cellType += value

        sample.antibody = antibody
        for value in sample.antibody.values():
            if re.search("h3k4me3", value, flags=re.IGNORECASE):
                ab_found = True
                break

        sample.treatment = treatment
        sample.disease = disease
        sample.cellLine = cellLine
        sample.genotype = genoType
        sample.tissue = tissue
        sample.title_found = title_found
        sample.ab_found = ab_found
        if title_found or ab_found:
            sample.title_ab = True

        if sample.organism == "Homo sapiens" and (sample.SRA != None and sample.SRA.strip() != "") and \
                sample.InstrumentID.startswith('Illu') and sample.libraryStrategy.lower() == "chip-seq":
            if sample.title_ab == True:
                # AllUniqueGSEs = AllUniqueGSEs.union(sample.series)
                Human_Samples[sample.id] = sample
                # HumanWithH3K4me3Download.add(sample.SRA)
            relatedSamples[sample.id] = sample
            for gse in sample.series:
                groupByGSE[gse].add(sample.id)
        file.close()

        # for char in characteristics.keys():
        #     totalCharacteristicsName[char]+=1
        # if len(feature) != 0:
        #     samples[sampleName] = sample
        #     totalOrganismsName[sampleOrganism]+=1
        # else:
        #     notFeature[sampleName] = sample

    # with open("./GEOSearchXMLs/geoSoftParserResult", "w") as file:
    #     for value in samples.values():
    #         json.dump(value.__dict__, file)
    # with open("./GEOSearchXMLs/geoSoftParserNoFeatureResult", "w") as file:
    #     for value in notFeature.values():
    #         json.dump(value.__dict__, file)
    # #
    # with open("./GEOSearchXMLs/geoorganimsWithH3K4me3.csv", "wb") as csv_file:
    #     writer = csv.writer(csv_file)
    #     for key, value in totalOrganismsName.items():
    #         writer.writerow([key, value])
    # #
    # with open("./GEOSearchXMLs/geocharacteristics.csv", "wb") as csv_file:
    #     writer = csv.writer(csv_file)
    #     for key, value in totalCharacteristicsName.items():
    #         writer.writerow([key, value])
    #
    # with open("./GEOSearchXMLs/geocontainingMessageWithH3K4me3.csv", "wb") as csv_file:
    #     writer = csv.writer(csv_file)
    #     for key, value in featureMeasage.items():
    #         writer.writerow([key, value])
    #         writer.writerow(["     "])

    # with open("./GEOSearchXMLs/geoH3K4me3GSMList.csv", "wb") as csv_file:
    #     writer = csv.writer(csv_file)
    #     writer.writerow(
    #         ['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
    #          "H3K4me3_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab",
    #          "Feature in Title or Ab"])
    #     for sample in samples.values():
    #         writer.writerow(
    #             [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
    #              sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine, sample.cellType,
    #              sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found,
    #              sample.title_ab])

    csv_file = open("./HumanH3K4me3GSMList.csv", "wb")
    writer = csv.writer(csv_file)
    writer.writerow(
        ['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
         "H3K4me3_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title",
         "Feature in Ab", "Feature in Title or Ab"])
    for sample in Human_Samples.values():
        writer.writerow(
            [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
             sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine, sample.cellType,
             sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found,
             sample.title_ab])
    csv_file.close()


    csv_file =  open("./HumanH3K4me3RelatedSamples.csv", "wb")
    writer = csv.writer(csv_file)
    writer.writerow(
        ['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID",
         "Library Strategy",
         "H3K4me3_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype",
         "Antibody", "Feature in Title",
         "Feature in Ab", "Feature in Title or Ab"])
    for sample in relatedSamples.values():
        writer.writerow(
            [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
             sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine,
             sample.cellType,
             sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found,
             sample.ab_found,
             sample.title_ab])
    csv_file.close()

    csv_file = open("./allHumanH3K4Samples.csv", "wb")
    writer = csv.writer(csv_file)
    writer.writerow(
        ['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
         "H3K4me3_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title",
         "Feature in Ab", "Feature in Title or Ab"])
    for sample in Human_Samples.values():
        writer.writerow(
            [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
             sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine, sample.cellType,
             sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found,
             sample.title_ab])
    for sample in relatedSamples.values():
        writer.writerow(
            [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
             sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine,
             sample.cellType,
             sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found,
             sample.ab_found,
             sample.title_ab])
    csv_file.close()

    # with open("./GEOSearchXMLs/notFeature.csv", "w") as csv_file:
    #     writer = csv.writer(csv_file)
    #     writer.writerow(['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
    #          "H3K4me3_description", "Tissue", "Cell Line", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab"])
    #     for sample in notFeature.values():
    #         writer.writerow(
    #             [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
    #              sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine,
    #              sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found])

    # with open("./allUniqueGSEsHumanWithH3K4me3.txt", "w") as file:
    #     for value in AllUniqueGSEs:
    #         file.write(value+"\n")
    #
    # with open("./HumanWithH3K4me3Download.txt", "w") as file:
    #     for value in HumanWithH3K4me3Download:
    #         file.write(value+"\n")
    return groupByGSE, Human_Samples, relatedSamples, encodeGSE


def spliterFinder(title, keyword):
    # find the spliter in the title, and return the keywords index and the spliter
    spliter = None
    index = None
    title = title.lower()

    space = title.count(" ")
    underscore = title.count("_")
    if space > 0 or underscore >0:
        if space > underscore:
            spliter = " "
        else:
            spliter = "_"
    if spliter is None:
        if title.find(keyword.lower()) != -1:
            return None, 0
        else:
            return None, None
    elements = title.split(spliter)
    for i in range(len(elements)):
        if elements[i] == keyword.lower():
            index = i
            break
    return spliter, index


def Similarity(title1, keyword1, title2, keyword2):
    title1 = title1.lower().replace(keyword1.lower(), "")
    title2 = title2.lower().replace(keyword2.lower(), "")

    return SequenceMatcher(None, title1, title2).ratio()




if __name__ == "__main__":
    groupbyGSE, HumanSamples, relatedSamples, encodeGSE = SOFTQuickRelated("./QuickXMLs", False)

    # initiate the map of sample to input
    SampleToInput = defaultdict(set)
    scoreBoard = {}
    #get all the sample with key word in title
    titleCandidates = set()
    for key, value in HumanSamples.items():
        if value.title_found == True:
            titleCandidates.add(key)
    not_found = 0
    # get their related samples
    for candidate in titleCandidates:
        sample = HumanSamples[candidate]

        sample_spliter, keyword_index = spliterFinder(sample.title, "H3K4me3")

        encode = False
        for gse in sample.series:
            if gse in encodeGSE:
                encode = True
                break
        if encode:
            targetGSEs = encodeGSE.union(set(sample.series))
        else:
            targetGSEs = set(sample.series)

        candidate = None
        bestMatchID = None
        bestSimilarity = float("-inf")
        for gse in targetGSEs:
            for relatedSample in groupbyGSE[gse]:
                input_spliter, input_index = spliterFinder(relatedSamples[relatedSample].title, "input")
                igg_spliter, igg_index = spliterFinder(relatedSamples[relatedSample].title, "IgG")
                wce_spliter, wce_index = spliterFinder(relatedSamples[relatedSample].title, "wce")
                control_spliter, control_index = spliterFinder(relatedSamples[relatedSample].title, "control")

                hasInput = [True if x == sample_spliter and y == keyword_index else False for x, y in
                            [(input_spliter, input_index),
                             (igg_spliter, igg_index),
                             (wce_spliter, wce_index),
                             (control_spliter, control_index)]]
                if hasInput[0] is True:
                    related_keyword = "input"
                elif hasInput[1] is True:
                    related_keyword = "IgG"
                elif hasInput[2] is True:
                    related_keyword = "WCE"
                elif hasInput[3] is True:
                    related_keyword = "control"

                if sample_spliter != None and keyword_index != None:
                    if any(hasInput):
                        score = Similarity(sample.title, "H3K4me3", relatedSamples[relatedSample].title, related_keyword)
                        if score > bestSimilarity:
                            bestSimilarity = score
                            bestMatchID = relatedSamples[relatedSample].id
                elif keyword_index == 0:
                    if relatedSamples[relatedSample].title.find("input") != -1:
                        SampleToInput[sample.id].add(relatedSamples[relatedSample].id)

                    elif relatedSamples[relatedSample].title.lower().find("wce") != -1:
                        SampleToInput[sample.id].add(relatedSamples[relatedSample].id)

                    elif relatedSamples[relatedSample].title.find("IgG") != -1:
                        SampleToInput[sample.id].add(relatedSamples[relatedSample].id)

                    elif relatedSamples[relatedSample].title.find("control") != -1:
                        SampleToInput[sample.id].add(relatedSamples[relatedSample].id)

                    # else:
                    #     for value in relatedSamples[relatedSample].characteristics.values():
                    #         if value.find("input"):
                    #             candidate = relatedSamples[relatedSample].id
                    #             break
                    #         elif value.lower().find(" wce "):
                    #             candidate = relatedSamples[relatedSample].id
                    #             break
                    #         elif value.lower().find("whole cell extract"):
                    #             candidate = relatedSamples[relatedSample].id
                    #             break


        if bestMatchID:
            SampleToInput[sample.id].add(bestMatchID)
	    scoreBoard[sample.id] = bestSimilarity

            #print bestSimilarity
	else:
	    not_found += 1

    for key, value in HumanSamples.items():
        if value.title_found != True:
            sample = value
            targetGSEs = set(sample.series)
            for gse in targetGSEs:
                for relatedSample in groupbyGSE[gse]:
                    for v in relatedSamples[relatedSample].characteristics.values():
                        if v.find("input"):
                            SampleToInput[sample.id].add(relatedSamples[relatedSample].id)
                            break
                        elif value.lower().find(" wce "):
                            SampleToInput[sample.id].add(relatedSamples[relatedSample].id)
                            break
                        elif value.lower().find("whole cell extract"):
                            SampleToInput[sample.id].add(relatedSamples[relatedSample].id)
                            break
                    if len(SampleToInput[sample.id]) > 0:
                        break
                if len(SampleToInput[sample.id]) > 0:
                    break
            print -1
            if len(SampleToInput[sample.id]) == 0:
                not_found+=1
    print not_found


    output = open("./H3K4me3_Sample_To_Input_BX_search.csv", "w")
    for key, value in SampleToInput.items():
        writer = csv.writer(output)
        row = [key]+[HumanSamples[key].title]+[scoreBoard[key]]
        # print value
        for id in value:
            row += [id]+[relatedSamples[id].title]
        writer.writerow(row)
    output.close()










