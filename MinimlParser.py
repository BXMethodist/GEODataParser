from xml.dom import minidom
from GSM import GSM
import os
import re
import json
from collections import defaultdict
import csv

def MinimlXMLParser(cwd=None):
    if cwd == None:
        cwd = os.getcwd()+'/xml'

    samples = {}
    featureOrganismsName = defaultdict(int)
    featureCharacteristicsName = defaultdict(int)
    featureMeasage = defaultdict(list)

    # cwd = os.getcwd()   #for test only
    file = open("xmlParserResult", "w")
    n = 0


    for filename in os.listdir(cwd):

        # filename = "GSE82260_family.xml"    #for test only
        # if not re.search("GSE[0-9]+", filename):
        #     continue
        # print filename



        doc = minidom.parse(cwd+"/"+filename)

        platFormName = doc.getElementsByTagName("Platform")[0].getAttribute("iid")

        sampleList = doc.getElementsByTagName("Sample")

        for sampleNode in sampleList:
            h3k4me3 = False
            # use nodeType to get the tpyes of node. See the constructor of xml.dom's Node class. Here the result is ElementNode.
            sampleName = sampleNode.getAttribute("iid")
            sample = GSM(sampleName)

            sample.platForm = platFormName

            channel = sampleNode.getElementsByTagName('Channel')
            characteristics = channel[0].getElementsByTagName('Characteristics')

            if len(characteristics) <= 1:
                sample.singleCharacteristics = None
            else:
                i = 0
                for char in characteristics:
                    # print type(char)
                    try:
                        charName = char._attrs.values()[0].nodeValue.strip()
                    except:
                        charName = "characteristic"+str(i)
                        i+=1
                    charValue = char.childNodes[0].nodeValue.strip()
                    if re.search('h3k4me3', charValue, flags=re.IGNORECASE):
                        h3k4me3 = True
                        featureMeasage[charName].append(sampleName + "   " + charValue)
                    sample.characteristics[charName] = charValue

            organism = channel[0].getElementsByTagName('Organism')
            organismName = organism[0].childNodes[0].nodeValue.strip()
            sample.organism = organismName

            try:
                libraryStrategy = sampleNode.getElementsByTagName("Library-Strategy")[0]
                sample.libraryStrategy = libraryStrategy.childNodes[0].nodeValue.lower().strip()
            except:
                pass

            title = sampleNode.getElementsByTagName("Title")[0]
            titleName = title.childNodes[0].nodeValue.strip()

            if re.search("h3k4me3", titleName, flags=re.IGNORECASE):
                h3k4me3 = True
                featureMeasage['Title'].append(titleName)
            sample.title = titleName

            # librarySelection = doc.getElementsByTagName("Library-Selection")[0]
            # sample.librarySelection = librarySelection.childNodes[0].nodeValue.lower().strip()

            supplementaryData = sampleNode.getElementsByTagName("Supplementary-Data")
            for data in supplementaryData:
                dataName = data._attrs.values()[0].nodeValue.strip()
                dataValue = data.childNodes[0].nodeValue.strip()
                sample.supplementaryData[dataName] = dataValue

            relations = sampleNode.getElementsByTagName("Relation")
            for relation in relations:
                if relation.getAttribute("type") == "SRA" :
                    sample.SRAurl = relation.getAttribute("target")
                    sample.SRA = sample.SRAurl[sample.SRAurl.find("=")+1:]
            if h3k4me3:
                if sampleName not in samples:
                    featureOrganismsName[sample.organism.lower()] += 1
                    for charName in sample.characteristics.iterkeys():
                        featureCharacteristicsName[charName] += 1
                    n += 1
                    samples[sampleName] = sample
                else:
                    # print sampleName, GSEName, samples[sampleName].series
                    pass

        # if samples.values() != []:
        #     break
    # print n
    # featureOrganismsName[samples['GSM1003803'].organism.lower()] -= 1
    # featureOrganismsName[samples['GSM1003805'].organism.lower()] -= 1
    # del samples['GSM1003803']
    # del samples['GSM1003805']
    with open("addxmlParserResult", "a") as file:
        for value in samples.values():
            # print sample.SRA
            json.dump(value.__dict__, file)

    with open("addorganimsWithH3K4me3.csv", "wb") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in featureOrganismsName.items():
            writer.writerow([key, value])

    with open("addcharacteristicsWithH3K4me3.csv", "wb") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in featureCharacteristicsName.items():
            writer.writerow([key, value])

    with open("addcontainingMessageWithH3K4me3.csv", "wb") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in featureMeasage.items():
            writer.writerow([key, value])
            writer.writerow(["     "])

    with open("addH3K4me3GSMList.csv", "wb") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['GSM_ID', "Organism", "GPL_ID", "SRA_ID", "Library Strategy"])
        for sample in samples.values():
            writer.writerow([sample.id, sample.organism, sample.series, sample.platForm, sample.SRA, sample.libraryStrategy])
    return samples, featureOrganismsName, featureCharacteristicsName

# for key, value in MinimlXMLParser().iteritems():
#     print value.title
samples, organismsName, characteristicsName = MinimlXMLParser("./GSMXMLs")
# values.sort(key=lambda x:x.title)




print "There are total ", len(samples), " samples"
print "From ", len(organismsName), " different organisms"
print organismsName
print "Has ", len(organismsName), " different Characteristics"
print characteristicsName