from xml.dom import minidom
from GSM import GSM
import os
import re
import json
from collections import defaultdict
import csv

def SOFTQuickParser(cwd=None, geo=True):
    if cwd == None:
        return

    # print len(map)
    samples = {}

    Human_Samples = {}

    totalOrganismsName = defaultdict(int)
    # totalCharacteristicsName = defaultdict(int)
    # featureMeasage = defaultdict(list)

    notFeature = {}


    geoGSMs = set()
    file =  open("uniqueGSM_GEOsearch.txt", "r")
    for line in file.readlines():
        geoGSMs.add(line.strip())
    file.close()

    for filename in os.listdir(cwd):
        if not filename.startswith("GSM"):
            continue

        if geo and filename[:-4] not in geoGSMs:
            continue

        sampleName = filename[:-4]
        sampleTitle = ""
        sampleType = ""
        sampleLibraryStrategy=''
        sampleOrganism=''
        samplePlatForm=''
        sampleInstrumentID=''

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
            if line.startswith("!Sample_title"):
                sampleTitle = line[line.find("=")+1:].strip()
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
                Human_Samples[sample.id] = sample
        file.close()

        # for char in characteristics.keys():
        #     totalCharacteristicsName[char]+=1
        if len(feature) != 0:
            samples[sampleName] = sample
            totalOrganismsName[sampleOrganism]+=1
        else:
            notFeature[sampleName] = sample

    # print totalOrganismsName

    # with open("./GEOSearchXMLs/geoSoftParserResult", "w") as file:
    #     for value in samples.values():
    #         json.dump(value.__dict__, file)
    # with open("./GEOSearchXMLs/geoSoftParserNoFeatureResult", "w") as file:
    #     for value in notFeature.values():
    #         json.dump(value.__dict__, file)
    # #
    if geo:
        outputOrganism = "./"+"GEOsearch"+"organimsWithH3K4me3.csv"
        outputHuman = "./"+"GEOsearch"+"humanWithH3K4me3.csv"
        outputSample = "./"+"GEOsearch"+"sampleWithH3K4me3.csv"
        outputNoFeature = "./"+"GEOsearch"+"noWithH3K4me3.csv"

    else:
        outputOrganism = "./" + "organimsWithH3K4me3.csv"
        outputHuman = "./" + "humanWithH3K4me3.csv"
        outputSample = "./" + "sampleWithH3K4me3.csv"
        outputNoFeature = "./" + "noWithH3K4me3.csv"


    with open(outputOrganism, "wb") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in totalOrganismsName.items():
            writer.writerow([key, value])

    with open(outputSample, "wb") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            ['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
             "H3K4me3_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab",
             "Feature in Title or Ab"])
        for sample in samples.values():
            writer.writerow(
                [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
                 sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine, sample.cellType,
                 sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found,
                 sample.title_ab])

    csv_file = open(outputHuman, "wb")
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


    # csv_file =  open("./HumanH3K4me3RelatedSamples.csv", "wb")
    # writer = csv.writer(csv_file)
    # writer.writerow(
    #     ['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID",
    #      "Library Strategy",
    #      "H3K4me3_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype",
    #      "Antibody", "Feature in Title",
    #      "Feature in Ab", "Feature in Title or Ab"])
    # for sample in relatedSamples.values():
    #     writer.writerow(
    #         [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
    #          sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine,
    #          sample.cellType,
    #          sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found,
    #          sample.ab_found,
    #          sample.title_ab])
    # csv_file.close()

    with open(outputNoFeature, "wb") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
             "H3K4me3_description", "Tissue", "Cell Line", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab"])
        for sample in notFeature.values():
            writer.writerow(
                [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
                 sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine,
                 sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found])

    # with open("./allUniqueGSEsHumanWithH3K4me3.txt", "w") as file:
    #     for value in AllUniqueGSEs:
    #         file.write(value+"\n")
    #
    # with open("./HumanWithH3K4me3Download.txt", "w") as file:
    #     for value in HumanWithH3K4me3Download:
    #         file.write(value+"\n")
    return

# for key, value in MinimlXMLParser().iteritems():
#     print value.title
# SOFTQuickParser("./GEOSearchXMLs")
# values.sort(key=lambda x:x.title)




# print "There are total ", len(samples), " samples"

# print "There are total", len(noFeature), " samples don't have features"
# print organismsName
# print "Has ", len(organismsName), " different Characteristics"
# print characteristicsName

SOFTQuickParser("/home/tmhbxx3/scratch/XMLhttp/QuickXMLs", geo=True)