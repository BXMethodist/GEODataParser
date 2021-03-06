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
import re, urllib2, gc, json, sqlite3, contextlib, pandas as pd
from collections import defaultdict
from difflib import SequenceMatcher
from update import GSE_info
from multiprocessing import Queue, Process
from sample import GSM


def get_WebInfo(url, count):
    with contextlib.closing(urllib2.urlopen(url)) as web:
        info = web.readlines()
    web.close()
    del web
    if count % 50 == 0:
        gc.collect()
    return info


def get_MetaInfo(db, sample, count):
    if db is not None:
        query = db.execute('select GSM_ID, MetaData from GSM where GSM_ID = "' + sample.id + '"').fetchall()

        if len(query) == 0:
            info = get_WebInfo(sample.url, count)
        else:
            try:
                info = json.loads(query[0][1])
            except:
                print "sample id is", sample.id
                print "This sample has problem !!!!!!", query[0][0]
                print "\t"
                print "\t"
                info = get_WebInfo(sample.url, count)
    else:
        info = get_WebInfo(sample.url, count)
    return info


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def SOFTQuickRelated(featured_samples, output_type, type_seq, GSEGSM_map, encode_remove, encodeGSE, cwd, process):
    relatedGSEs = []

    for key, value in featured_samples.iteritems():
        relatedGSEs += value.series

    relatedGSEs = set(relatedGSEs)
    if encode_remove:
        allrelatedGSEs = relatedGSEs
    else:
        allrelatedGSEs = encodeGSE.union(relatedGSEs)
    allrelatedGSEs = list(allrelatedGSEs)

    allrelatedGSMs = set()

    for gse in allrelatedGSEs:
        if gse in GSEGSM_map:
            allrelatedGSMs = allrelatedGSMs.union(GSEGSM_map[gse])
        else:
            allrelatedGSMs = allrelatedGSMs.union(GSE_info(gse)[0])

    allrelatedGSMs = list(allrelatedGSMs)

    # print "total ", len(allrelatedGSMs), " samples found for potential input"

    chunksize = len(allrelatedGSMs)/(process-1) if process > 1 else len(allrelatedGSMs)

    # print "chunksize is ", chunksize

    queue = Queue()
    processes = []
    for i in range(process):
        cur_relatedGSMs = allrelatedGSMs[i * chunksize:(i + 1) * chunksize]
        p = Process(target=related_sample_info, args=(cur_relatedGSMs, queue, output_type, type_seq, cwd, featured_samples))
        processes.append(p)
        p.start()

    relatedSamples = {}
    groupByGSE = defaultdict(set)

    for i in range(process):
        cur_relatedSamples, cur_groupByGSE = queue.get()
        relatedSamples.update(cur_relatedSamples)
        for key, value in cur_groupByGSE.items():
            groupByGSE[key]=groupByGSE[key].union(value)
    for p in processes:
        p.join()
    return groupByGSE, encodeGSE, relatedSamples


def related_sample_info(cur_relatedGSMs, queue, output_type, type_seq, cwd, featured_samples):
    # print "Process id is ", os.getpid()

    relatedSamples = {}
    groupByGSE = defaultdict(set)

    count = 0

    if cwd is None:
        db = None
    else:
        db = sqlite3.connect(cwd)
        # db.text_factory = str

    for filegsm in cur_relatedGSMs:
        characteristics = defaultdict(str)
        supplementaryData = defaultdict(str)
        relations = defaultdict(str)
        sampleSeriesID = set()

        antibody = {}
        title_found = False
        ab_found = False

        sampleTitle = ""
        sampleType = ""
        sampleLibraryStrategy = ""
        sampleOrganism = ""
        samplePlatForm = ""
        sampleInstrumentID= ""

        sample = GSM(filegsm)

        info = get_MetaInfo(db, sample, count)

        count +=1

        if info is None:
            continue

        for line in info:
            if line.startswith("!Sample_title"):
                sampleTitle = line[line.find("=")+1:].strip()
                if sampleTitle.find(";") != -1:
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

        cellLine = ""
        cellType = ""
        tissue = ""

        for key, value in characteristics.items():
            if key.lower() in ["chip antibody", "chip", "antigen", "antibody", "antibodies", "chip antibodies",
                               "antibody name", "antibody target", "target", "antibody/capture", "antibody/vendor/catalog#",
                               "chip ab", "chip antibody1", "chip antibody2", "chip-antibody", "chip_antibody",
                               "chip-seq antibody", "chip_antibodies", "chip-antibodies", "histone mark", "epigenetic feature",
                               "histone modification", "antibody antibodydescription", "chip antibody (epitope/name)",
                               "factor", "chip antibody/mbd affinity column", "chip/dip antibody", "antibody epiptope",
                               "antibody source", 'modification', "antibody (vendor': ' catalog#, or reference)",
                                "purification antibody", "antibody/details", "antibody epiptope",
                               "antibody information", "chip antibody / digestive enzyme", "chip antiboy",
                               "ip antibody", "chip antibody target", "modification", "histone", "enrichment procedure",
                               "antibody (vendor': ' catalog#, or reference)", "developmental stage/condition/extract protocol",
                               "antibody source"] \
                    or re.search('antibody epitope|using[\w\s]+antibod|immunoprecipitat', key, flags=re.IGNORECASE):
                if key in antibody:
                    antibody[key] += ", " + value
                else:
                    antibody[key] = value

            if key.lower() in ["cell line",  "cell", "cells pointed by barcodes", "cell line background",
                           "chicken line", "line"]:
                cellLine += value

            if key.lower() in ["cell_type", "cell-type", "cell type", "cell lineage"]:
                cellType += value

            if key.lower() in ["tissue", "body part", "body site", "tissue type"]:
                tissue += value

        sample.antibody = antibody
        sample.cellLine = cellLine
        sample.cellType = cellType
        sample.tissue = tissue

        sample.title_found = title_found
        sample.ab_found = ab_found
        if title_found or ab_found:
            sample.title_ab = True

        if (sample.organism == output_type or output_type is None) and (sample.SRA != None or sample.SRA.strip() != "") and \
                sample.InstrumentID.startswith('Illu') and (sample.libraryStrategy.lower() == type_seq or type_seq is None)\
                and sample.id not in featured_samples:
            relatedSamples[sample.id] = sample
            for gse in sample.series:
                groupByGSE[gse].add(sample.id)
    queue.put((relatedSamples, groupByGSE))
    if db is not None:
        db.close()
    return


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


def Similarity(title1, keyword1, title2, keyword2, ignorecase):
    score = SequenceMatcher(None, title1, title2).ratio()

    if ignorecase:
        title1 = title1.lower().replace(keyword1.lower(), "").replace("chip-seq", "")
        title2 = title2.lower().replace(keyword2.lower(), "").replace("chip-seq", "")

        title1 = re.sub(r'rep[0-9]*', '', title1)
        title2 = re.sub(r'rep[0-9]*', '', title2)
    else:
        title1 = title1.replace(keyword1, "").lower().replace("chip-seq", "")
        title2 = title2.replace(keyword2, "").lower().replace("chip-seq", "")

        title1 = re.sub(r'rep[0-9]*', '', title1)
        title2 = re.sub(r'rep[0-9]*', '', title2)

    score_replace = SequenceMatcher(None, title1, title2).ratio()

    if score > score_replace:
        return score, score_replace
    else:
        return score_replace, score
    #return score_replace

def keyword(message, features, features_begin, ignorecase):
    if ignorecase:
        for feature in features:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.search(feature, message, flags=re.IGNORECASE):
                return feature
        for feature in features_begin:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.match(feature, message, flags=re.IGNORECASE):
                return feature
    else:
        for feature in features:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.search(feature, message):
                return feature
        for feature in features_begin:
            feature = feature.replace(" ","")
            feature = feature.replace("-", "")
            feature = feature.replace("_", "")
            if re.match(feature, message):
                return feature


def Character_Similarity(sample1, sample2):
    title1 = re.sub("\d", "", sample1.title)
    title2 = re.sub("\d", "", sample2.title)
    score = similar(title1, title2)
    for key, value in sample1.characteristics.items():
        if key in sample2.characteristics:
            score += similar(value, sample2.characteristics[key])

    # print sample1.id, sample2.id, score
    return score


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


def has_antibody(sample, keyword):
    for v in sample.antibody.values():
        if v.find(keyword) != -1:
            return True
    return False


def equal_antibody(sample, keyword):
    for v in sample.antibody.values():
        if v == keyword:
            return True
    return False


def isInput(sample, feature_key_word):
    non_capital_keywords = ['input','wce']

    if feature_key_word.find("H3K") == -1:
        capital_keywords = ['IgG']
    else:
        capital_keywords = ['IgG', '_H3_', " H3"]

    for c in non_capital_keywords:
        if sample.title.lower().find(c) != -1:
            return True, c
        if has_antibody(sample, c):
            return True, c

    for n in capital_keywords:
        if n ==' H3' and sample.title[-3:] == ' H3':
            return True, 'H3'
        elif sample.title.find(n) != -1 and n !=' H3':
            if n == "_H3_":
                return True, 'H3'
            return True, n
    return False, ""


def input_finder(output_surffix, output_path, HumanSamples, groupByGSE, encodeGSE, relatedSamples,
                 features, features_begin, ignorecase, output_type):
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

    # print "title and none title", len(titleCandidates), len(noneTitle)
    not_found = 0
    # get their related samples
    for candidate in titleCandidates:
        sample = HumanSamples[candidate]

        feature_key_word = keyword(sample.title, features, features_begin, ignorecase)

        # print feature_key_word

        encode = False
        for gse in sample.series:
            if gse in encodeGSE:
                encode = True
                break
        if encode:
            targetGSEs = encodeGSE.union(set(sample.series))
        else:
            targetGSEs = set(sample.series)

        bestMatchID = set()
        bestSimilarity = float("-inf")
        bestSimilarity2 = float("-inf")
        input_keyword = ""

        for gse in targetGSEs:
            for relatedSample in groupByGSE[gse]:
                score = None
                boo, word = isInput(relatedSamples[relatedSample], feature_key_word)
                # print relatedSamples[relatedSample].title, boo, word
                if boo \
                        and sample.cellLine == relatedSamples[relatedSample].cellLine \
                        and sample.cellType == relatedSamples[relatedSample].cellType \
                        and sample.tissue == relatedSamples[relatedSample].tissue:

                    score1, score2 = Similarity(sample.title, feature_key_word, relatedSamples[relatedSample].title,
                                       word, ignorecase)

                    if score1 > bestSimilarity:
                        bestSimilarity = score1
                        bestSimilarity2 = score2
                        bestMatchID = set()
                        bestMatchID.add(relatedSamples[relatedSample].id)
                        input_keyword = word
                    elif score1 == bestSimilarity:
                        if input_keyword == 'H3' and word != 'H3':
                            bestMatchID = set()
                            bestMatchID.add(relatedSamples[relatedSample].id)
                            input_keyword = word
                            bestSimilarity = score1
                            bestSimilarity2 = score2
                        elif input_keyword != 'H3' and word == 'H3':
                            pass
                        else:
                            if score2 > bestSimilarity2:
                                bestSimilarity2 = score2
                                bestMatchID = set()
                                bestMatchID.add(relatedSamples[relatedSample].id)

        if bestMatchID:
            FirstSampleToInput[sample.id] = FirstSampleToInput[sample.id].union(bestMatchID)
        else:
            not_found += 1

    for key in noneTitle:
        sample = HumanSamples[key]
        targetGSEs = set(sample.series)

        best_char_score = 0
        best_id = set()

        for gse in targetGSEs:
            for relatedSample in groupByGSE[gse]:
                char_score = None
                for v in relatedSamples[relatedSample].antibody.values():
                    if v.find("input")!= -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                    elif v.lower().find(" wce ") != -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                    elif v.lower().find("whole cell extract") != -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                    elif v.find("IgG") != -1 and sample.id != relatedSamples[relatedSample].id \
                            and sample.cellLine == relatedSamples[relatedSample].cellLine:
                        char_score = Character_Similarity(sample, relatedSamples[relatedSample])
                if char_score > best_char_score:
                    best_id = set()
                    best_id.add(relatedSamples[relatedSample].id)
                    best_char_score = char_score
                elif char_score == best_char_score:
                    best_id.add(relatedSamples[relatedSample].id)

        if best_id:
            ThirdSampleToInput[sample.id] = best_id
        else:
            not_found+=1

    # print not_found
    output_type = output_type.replace(" ", "_")
    output1 = output_path + "First_" + output_surffix + "_" + output_type + "_Sample_To_Input.csv"
    output3 = output_path + "Third_" + output_surffix + "_" + output_type +"_Sample_To_Input.csv"

    table = []
    for key, value in FirstSampleToInput.items():
        row = [key] + [HumanSamples[key].title] + [HumanSamples[key].antibody]
        for id in value:
            row += [id] + [relatedSamples[id].title] + [relatedSamples[id].antibody]
        table.append(row)

    df = pd.DataFrame(table, columns=None)
    # df.to_csv(output1, sep=',', encoding='utf-8')

    table = []

    for key, value in ThirdSampleToInput.items():
        row = [key] + [HumanSamples[key].antibody]
        for id in value:
            row += [id] + [relatedSamples[id].antibody]
        table.append(row)
    df = pd.DataFrame(table, columns=None)
    # df.to_csv(output3, sep=',', encoding='utf-8')

    return FirstSampleToInput, ThirdSampleToInput



