#### TO DO, create new database version


from GSM import GSM
import os
import re
from collections import defaultdict
import csv
from difflib import SequenceMatcher
import sqlite3

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def SOFTQuickRelated(featured_samples, cwd):
    if cwd == None:
        return

    relatedSamples = {}

    groupByGSE = defaultdict(set)


    db = sqlite3.connect('/home/tmhbxx3/archive/GEO_MetaDatabase/geoMetaData.db')
    db.text_factory = str

    relatedGSEs = []

    featureGSMs = set()

    for key, value in featured_samples:
        featureGSMs.add(key)
        relatedGSEs += value.series

    relatedGSEs = set(relatedGSEs)

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
        print n
        n += 998
        if n < len(allrelatedGSEs):
            block = allrelatedGSEs[n-998: n]
        else:
            block = allrelatedGSEs[n-998: -1]
            keep = False
        query = db.execute("SELECT distinct GSM_ID FROM GSEtoGSM WHERE GSE_ID IN (" + ",".join("?" * len(block)) + ")", block).fetchall()

        for gsmid in query:
            if gsmid not in featureGSMs:
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

        sample.title_found = title_found
        sample.ab_found = ab_found
        if title_found or ab_found:
            sample.title_ab = True

        if sample.organism == "Homo sapiens" and (sample.SRA != None or sample.SRA.strip() != "") and \
                sample.InstrumentID.startswith('Illu') and sample.libraryStrategy.lower() == "chip-seq":
            relatedSamples[sample.id] = sample
            for gse in sample.series:
                groupByGSE[gse].add(sample.id)
        file.close()

    return groupByGSE, encodeGSE, relatedSamples


def spliterFinder(title, keyword):
    # find the spliter in the title, and return the keywords index and the spliter
    spliter = None
    index = None
    title = title.lower()

    space = title.count(" ")
    underscore = title.count("_")
    hyphen = title.count("-")
    choices = [" ", "_", "-"]
    counts = [space, underscore, hyphen]
    if space > 0 or underscore >0 or hyphen >0:
        max_count = 0
        choice = None
        for i in range(len(choices)):
            if counts[i] > max_count:
                max_count = counts[i]
                choice = i
        spliter = choices[choice]

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

def keyword(message, features, features_begin, ignorecase):
    if ignorecase:
        for feature in features:
            if re.search(feature, message, flags=re.IGNORECASE):
                return feature
        for feature in features_begin:
            if re.match(feature, message, flags=re.IGNORECASE):
                return feature
    else:
        for feature in features:
            if re.search(feature, message):
                return feature
        for feature in features_begin:
            if re.match(feature, message):
                return feature

def input_finder(output_surffix, HumanSamples, groupByGSE, encodeGSE, relatedSamples,
                 features, features_begin, ignorecase):
    FirstSampleToInput = defaultdict(set)

    ThirdSampleToInput = defaultdict(set)

    # get all the sample with key word in title
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

        feature_key_word = keyword(sample.title, features, features_begin, ignorecase)

        sample_spliter, keyword_index = spliterFinder(sample.title, feature_key_word)

        encode = False
        for gse in sample.series:
            if gse in encodeGSE:
                encode = True
                break
        if encode:
            targetGSEs = encodeGSE.union(set(sample.series))
        else:
            targetGSEs = set(sample.series)

        bestMatchID = None
        bestSimilarity = float("-inf")

        related_keyword = None

        for gse in targetGSEs:
            for relatedSample in groupByGSE[gse]:
                score = None
                if keyword_index != None:
                    if relatedSamples[relatedSample].title.lower().find("input") != -1:
                        score = Similarity(sample.title, feature_key_word, relatedSamples[relatedSample].title,
                                           "input")
                        if related_keyword == None or related_keyword == "input":
                            related_keyword = "input"
                            if score > bestSimilarity:
                                bestSimilarity = score
                                bestMatchID = relatedSamples[relatedSample].id
                        elif related_keyword == "wce":
                            if score > bestSimilarity:
                                bestSimilarity = score
                                bestMatchID = relatedSamples[relatedSample].id
                                related_keyword = "input"
                        else:
                            related_keyword = "input"
                            bestSimilarity = score
                            bestMatchID = relatedSamples[relatedSample].id

                    elif relatedSamples[relatedSample].title.lower().find("wce") != -1:
                        score = Similarity(sample.title, feature_key_word, relatedSamples[relatedSample].title,
                                           "wce")
                        if related_keyword == None or related_keyword == "wce":
                            related_keyword = "wce"
                            if score > bestSimilarity:
                                bestSimilarity = score
                                bestMatchID = relatedSamples[relatedSample].id
                        elif related_keyword == "input":
                            if score > bestSimilarity:
                                bestSimilarity = score
                                bestMatchID = relatedSamples[relatedSample].id
                                related_keyword = "wce"
                        elif related_keyword == "IgG" or related_keyword == "control":
                            related_keyword = "wce"
                            bestSimilarity = score
                            bestMatchID = relatedSamples[relatedSample].id
                    elif relatedSamples[relatedSample].title.lower().find("IgG") != -1:
                        score = Similarity(sample.title, feature_key_word, relatedSamples[relatedSample].title,
                                           "IgG")
                        if related_keyword == None or related_keyword == "IgG":
                            related_keyword = "IgG"
                            if score > bestSimilarity:
                                bestSimilarity = score
                                bestMatchID = relatedSamples[relatedSample].id
                        elif related_keyword == "control":
                            related_keyword = "IgG"
                            bestSimilarity = score
                            bestMatchID = relatedSamples[relatedSample].id
                    elif relatedSamples[relatedSample].title.lower().find("control") != -1:
                        score = Similarity(sample.title, feature_key_word, relatedSamples[relatedSample].title,
                                           "control")
                        if related_keyword == None or related_keyword == "control":
                            related_keyword = "control"
                            if score > bestSimilarity:
                                bestSimilarity = score
                                bestMatchID = relatedSamples[relatedSample].id

        if bestMatchID:
            FirstSampleToInput[sample.id].add(bestMatchID)
        else:
            not_found += 1

    for key in noneTitle:
        sample = HumanSamples[key]
        targetGSEs = set(sample.series)
        for gse in targetGSEs:
            for relatedSample in groupByGSE[gse]:
                for v in relatedSamples[relatedSample].characteristics.values():
                    if v.find("input") and sample.id != relatedSamples[relatedSample].id \
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
            not_found += 1

    print not_found

    output1 = "./First_" + output_surffix + "_Sample_To_Input.csv"
    output3 = "./Third_" + output_surffix + "_Sample_To_Input.csv"

    output = open(output1, "w")
    for key, value in FirstSampleToInput.items():
        writer = csv.writer(output)
        row = [key] + [HumanSamples[key].title]
        # print value
        for id in value:
            row += [id] + [relatedSamples[id].title]
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

    return FirstSampleToInput, ThirdSampleToInput

def has_features(message, features, features_begin, ignorecase):
    if ignorecase:
        for feature in features:
            if re.search(feature, message, flags=re.IGNORECASE):
                return True
        for feature in features_begin:
            if re.match(feature, message, flags=re.IGNORECASE):
                return True
    else:
        for feature in features:
            if re.search(feature, message):
                return True
        for feature in features_begin:
            if re.match(feature, message):
                return True
    return False
