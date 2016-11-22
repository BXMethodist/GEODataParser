from xml.dom import minidom
from GSM import GSM
from xml.parsers import expat
import os
import re
import json
from collections import defaultdict
import csv

def GSEQuickXML(cwd=None):
    if cwd == None:
        cwd = os.getcwd()+'/XMLs'

    samples = []
    GSMs = set()
    wrongFormat = set()
    for filename in os.listdir("./XMLs/"):
        GSEName = filename[:filename.find(".xml")]
        try:
            doc = minidom.parse(cwd+"/"+filename, parser=expat.ParserCreate('UTF-8'))

            sampleList = doc.getElementsByTagName("Sample")

            for sampleNode in sampleList:
                sampleName = sampleNode.getAttribute("iid")
                sample = GSM(sampleName, GSEName)
                GSMs.add(sampleName)
                samples.append(sample)
        except:
            with open(cwd+"/"+filename) as file:
                for line in file.readlines():
                    if re.search('GSM[0-9]+', line):
                        for gsm in re.findall('GSM[0-9]+', line):
                            wrongFormat.add(gsm)



    with open("GSEQuickResult.csv", "wb") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['GSM_ID', "GSE_ID"])
        for sample in samples:
            writer.writerow([sample.id, sample.series])
    with open("ALLGSMsID.txt", "w") as file:
        for gsmId in GSMs:
            file.write(gsmId+"\n")
    with open("failedGSEXML.txt", "w") as file:
        for fail in wrongFormat:
            file.write(fail+"\n")

    return samples, GSMs

print GSEQuickXML()[1]