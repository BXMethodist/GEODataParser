from GSM import GSM
import os
import csv
from collections import defaultdict

def GSEQuickSOFT(cwd=None):
    if cwd == None:
        cwd = os.getcwd()+'/GSESoftQuick/'

    GSMtoGSE = defaultdict(set)
    GSEtoGSM = defaultdict(set)

    for filename in os.listdir(cwd):
        GSEName = filename[:filename.find(".txt")]

        with open(cwd+filename, "r") as file:
            for line in file.readlines():
                if line.startswith("!Series_sample_id"):
                    GSMid = line[line.find("GSM"):].strip()
                    GSMtoGSE[GSMid].add(GSEName)
                    GSEtoGSM[GSEName].add(GSMid)

    # with open("GSEQuickSOFTResult.csv", "wb") as csv_file:
    #     writer = csv.writer(csv_file)
    #     writer.writerow(['Sample_ID', "Series_ID"])
    #     for sample in samples.values():
    #         writer.writerow([sample.id, sample.series])
    GSMsWithFeauture = set()
    with open("allUniqueGSEsHumanWithH3K4me3.txt", "r") as file:
        for line in file.readlines():
            if len(GSEtoGSM[line.strip()]) == 0:
                print "wrong", line
            GSMsWithFeauture = GSMsWithFeauture.union(GSEtoGSM[line.strip()])


    with open("allUniqueGSMsHumanWithH3K4me3.txt", "w") as file:
        for line in GSMsWithFeauture:
            file.write('%s\n' %line)


    with open("GSMtoGSE.txt", "w") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in GSMtoGSE.items():
            writer.writerow([key, list(value)])

    with open("GSEtoGSM.txt", "w") as csv_file:
        writer = csv.writer(csv_file)
        for key, value in GSEtoGSM.items():
            writer.writerow([key, list(value)])


    return

GSEQuickSOFT()