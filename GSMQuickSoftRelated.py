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

    featureGSMs = list(featureGSMs)
    m = 0
    keep = True
    while keep:
        m += 998
        if m < len(featureGSMs):
            block = featureGSMs[m-998: m]
        else:
            block = featureGSMs[m-998: -1]
            keep = False
        query = db.execute("SELECT distinct GSE_ID FROM GSEtoGSM WHERE GSM_ID IN (" + ",".join("?" * len(block)) + ")",
                           block).fetchall()

        for gseid in query:
            relatedGSEs.add(gseid[0])

    print len(relatedGSEs)

    encodeGSE = set()
    allGSEs = db.execute('select distinct GSE_ID from GSE where organization = "ENCODE DCC"').fetchall()
    for gse in allGSEs:
        encodeGSE.add(gse[0])

    print len(encodeGSE)

    allrelatedGSEs = encodeGSE.union(relatedGSEs)
    allrelatedGSEs = list(allrelatedGSEs)

    allrelatedGSMs = set()

    n = 0
    keep = True
    while keep:
        n += 998
        if n < len(allrelatedGSEs):
            block = allrelatedGSEs[n-998: n]
        else:
            block = allrelatedGSEs[n-998: -1]
            keep = False
        query = db.execute("SELECT distinct GSM_ID FROM GSEtoGSM WHERE GSE_ID IN (" + ",".join("?" * len(block)) + ")", block).fetchall()

        for gsmid in query:
            allrelatedGSMs.add(gsmid[0])

    db.close()

    print len(allrelatedGSMs)

    allrelatedGSMs = list(allrelatedGSMs)

    for filegsm in allrelatedGSMs:
        filename = filegsm+".xml"


        file = open(cwd+'/'+filename, "r")
        characteristics = defaultdict(str)
        supplementaryData = defaultdict(str)
        relations = defaultdict(str)
        sampleSeriesID = set()
        feature = {}

        antibody = {}
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

        sample.antibody = antibody
        for value in sample.antibody.values():
            if re.search("h3k4me3", value, flags=re.IGNORECASE):
                ab_found = True
                break

        sample.title_found = title_found
        sample.ab_found = ab_found
        if title_found or ab_found:
            sample.title_ab = True

        if sample.organism == "Homo sapiens" and (sample.SRA != None or sample.SRA.strip() != "") and \
                sample.InstrumentID.startswith('Illu') and sample.libraryStrategy.lower() == "chip-seq":
            if sample.title_ab == True:
                Human_Samples[sample.id] = sample

            relatedSamples[sample.id] = sample
            for gse in sample.series:
                groupByGSE[gse].add(sample.id)
        file.close()

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
            return None, -1
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
    groupbyGSE, HumanSamples, relatedSamples, encodeGSE = SOFTQuickRelated("/home/tmhbxx3/scratch/XMLhttp/QuickXMLs", False)

    geo = False

    print "parser done!"
    print "groupbyGSE size is ", len(groupbyGSE)
    print "HumanSamples size is ", len(HumanSamples)
    print "relatedSamples size is ", len(relatedSamples)

    # initiate the map of sample to input
    FirstSampleToInput = defaultdict(set)

    SecondSampleToInput = defaultdict(set)

    ThirdSampleToInput = defaultdict(set)

    #get all the sample with key word in title
    titleCandidates = set()

    noneTitle = set()
    for key, value in HumanSamples.items():
        if value.title_found == True:
            titleCandidates.add(key)
        else:
            noneTitle.add(key)

    print "title and none title", len(titleCandidates), len(noneTitle)
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
                if sample_spliter != None and keyword_index != None:
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

                    if any(hasInput):
                        score = Similarity(sample.title, "H3K4me3", relatedSamples[relatedSample].title, related_keyword)
                        if score > bestSimilarity:
                            bestSimilarity = score
                            bestMatchID = relatedSamples[relatedSample].id
                elif keyword_index == -1:
                    if relatedSamples[relatedSample].title.find("input") != -1:
                        score = Similarity(sample.title, "H3K4me3", relatedSamples[relatedSample].title,
                                           "input")

                    elif relatedSamples[relatedSample].title.lower().find("wce") != -1:
                        score = Similarity(sample.title, "H3K4me3", relatedSamples[relatedSample].title,
                                           "wce")

                    elif relatedSamples[relatedSample].title.find("IgG") != -1:
                        score = Similarity(sample.title, "H3K4me3", relatedSamples[relatedSample].title,
                                           "IgG")
                    if score > bestSimilarity:
                        bestSimilarity = score
                        bestMatchID = relatedSamples[relatedSample].id


        if bestMatchID and keyword_index != -1:
            FirstSampleToInput[sample.id].add(bestMatchID)
        elif bestMatchID and keyword_index == -1:
            SecondSampleToInput[sample.id].add(bestMatchID)
        else:
            not_found+=1


    for key in noneTitle:
        sample = HumanSamples[key]
        targetGSEs = set(sample.series)
        for gse in targetGSEs:
            for relatedSample in groupbyGSE[gse]:
                for v in relatedSamples[relatedSample].characteristics.values():
                    if v.find("input") and sample.id != relatedSamples[relatedSample].id\
                            and relatedSamples[relatedSample].title.lower().find("h3k") == -1:
                        ThirdSampleToInput[sample.id].add(relatedSamples[relatedSample].id)
                        break
                    elif v.lower().find(" wce ") and sample.id != relatedSamples[relatedSample].id \
                            and relatedSamples[relatedSample].title.lower().find("h3k") == -1:
                        ThirdSampleToInput[sample.id].add(relatedSamples[relatedSample].id)
                        break
                    elif v.lower().find("whole cell extract") and sample.id != relatedSamples[relatedSample].id \
                            and relatedSamples[relatedSample].title.lower().find("h3k") == -1:
                        ThirdSampleToInput[sample.id].add(relatedSamples[relatedSample].id)
                        break
                if sample.id in ThirdSampleToInput:
                    break
            if sample.id in ThirdSampleToInput:
                break

        if not sample.id in ThirdSampleToInput:
            not_found+=1

    print not_found


    if geo:
        output1 = "./First_H3K4me3_Sample_To_Input.csv"
        output2 = "./Second_H3K4me3_Sample_To_Input.csv"
        output3 = "./Third_H3K4me3_Sample_To_Input.csv"
    else:
        output1 = "./First_H3K4me3_Sample_To_Input_my_search.csv"
        output2 = "./Second_H3K4me3_Sample_To_Input_my_search.csv"
        output3 = "./Third_H3K4me3_Sample_To_Input_my_search.csv"

    output = open(output1, "w")
    for key, value in FirstSampleToInput.items():
        writer = csv.writer(output)
        row = [key]+[HumanSamples[key].title]
        # print value
        for id in value:
            row += [id]+[relatedSamples[id].title]
        writer.writerow(row)
    output.close()

    output = open(output2, "w")
    for key, value in SecondSampleToInput.items():
        writer = csv.writer(output)
        row = [key]+[HumanSamples[key].title]
        # print value
        for id in value:
            row += [id]+[relatedSamples[id].title]
        writer.writerow(row)
    output.close()

    output = open(output3, "w")
    for key, value in ThirdSampleToInput.items():
        writer = csv.writer(output)
        row = [key] + [HumanSamples[key].title]
        # print value
        for id in value:
            row += [id] + [relatedSamples[id].title]
        writer.writerow(row)
    output.close()








