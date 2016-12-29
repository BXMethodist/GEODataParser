from xml.dom import minidom
from GSM import GSM
import os
import re
import json
from collections import defaultdict
import csv

def has_features(message, features):
    for feature in features:
        feature = feature.lower()
        if re.search(feature, message, flags=re.IGNORECASE):
            return True
    return False

def SOFTQuickParser(output_surfix, cwd=None, geo=True, *features):
    if cwd == None:
        return

    # print len(map)
    samples = {}

    Human_Samples = {}

    totalOrganismsName = defaultdict(int)

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
        target_feature = {}

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
                if has_features(sampleTitle, features):
                    target_feature["Title"] = sampleTitle
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
                if has_features(value, features):
                    target_feature[key] = value
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
        sample.features = target_feature
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
            if has_features(value, features):
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
        if len(target_feature) != 0:
            samples[sampleName] = sample
            totalOrganismsName[sampleOrganism]+=1
        else:
            notFeature[sampleName] = sample

    if geo:
        outputOrganism = "./"+"GEOsearch"+"organimsWith" + output_surfix +".csv"
        outputHuman = "./"+"GEOsearch"+"humanWith" + output_surfix + ".csv"
        outputSample = "./"+"GEOsearch"+"sampleWith" + output_surfix + ".csv"
        outputNoFeature = "./"+"GEOsearch"+"noWith" + output_surfix+ ".csv"

    else:
        outputOrganism = "./" + "organimsWith" + output_surfix+ ".csv"
        outputHuman = "./" + "humanWith" + output_surfix + ".csv"
        outputSample = "./" + "sampleWith" + output_surfix + ".csv"
        outputNoFeature = "./" + "noWith" + output_surfix + ".csv"


    with open(outputOrganism, "wb") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in totalOrganismsName.items():
            writer.writerow([key, value])

    with open(outputSample, "wb") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(
            ['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
             output_surfix + "_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab",
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
         output_surfix+"_description", "Tissue", "Cell Line", "Cell Type", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title",
         "Feature in Ab", "Feature in Title or Ab"])
    for sample in Human_Samples.values():
        writer.writerow(
            [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
             sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine, sample.cellType,
             sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found,
             sample.title_ab])
    csv_file.close()

    if geo:
        with open(outputNoFeature, "wb") as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(['Sample_ID', "Title", "Organism", "Series_ID", "GPL_ID", "Instrument Model", "SRA_ID", "Library Strategy",
                output_surfix+"_description", "Tissue", "Cell Line", "Disease", "Treatment", "Genotype", "Antibody", "Feature in Title", "Feature in Ab"])
            for sample in notFeature.values():
                writer.writerow(
                    [sample.id, sample.title, sample.organism, sample.series, sample.platForm, sample.InstrumentID,
                     sample.SRA, sample.libraryStrategy, sample.features, sample.tissue, sample.cellLine,
                     sample.disease, sample.treatment, sample.genotype, sample.antibody, sample.title_found, sample.ab_found])
    return

SOFTQuickParser("H3K27me3", "/home/tmhbxx3/scratch/XMLhttp/QuickXMLs", False, "H3K27me3", "K27me3")
