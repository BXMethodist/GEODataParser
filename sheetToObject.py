# This module is used for convert the parse sheet results to GSM objects
from GSM import GSM
import pandas as pd
import numpy as np
import ast
from collections import defaultdict
import csv

def sheetToObjects(path):
    table = pd.read_csv(path, header=0)

    GSE_index = table.columns.get_loc("Series_ID")
    GSE_ID_strs = table.values[:, GSE_index]
    GSE_IDs = set()
    for id in GSE_ID_strs:
        ids = id[1: -1].split(",")
        for i in ids:
            i = i.replace("'", "")
            GSE_IDs.add(i.strip())


    SRA_index = table.columns.get_loc("SRA_ID")
    SRA_links = table.values[:, SRA_index]

    with open("related_SRA_links.csv", "w") as file:
        for link in SRA_links:
            file.write(link + "\n")
    with open("related_GSE_IDs.csv", "w") as file:
        for id in GSE_IDs:
            file.write(id + "\n")

    samples = {}

    for rowNumber in range(table.shape[0]-1):
        sample = GSM(table.ix[rowNumber, "Sample_ID"])
        sample.title = table.ix[rowNumber, "Title"]
        sample.organism = table.ix[rowNumber, "Organism"]
        sample.series = ast.literal_eval(table.ix[rowNumber, "Series_ID"])
        sample.cellLine = table.ix[rowNumber, "Cell Line"]
        sample.cellType = table.ix[rowNumber, "Cell Type"]

        samples[sample.id] = sample

    return samples



def getRelatedGSMs(path="./related_GSE_IDs.csv", map="./GSEtoGSM.txt"):
    GSE = open(path, "r")
    GSEs = [line.rstrip() for line in GSE.readlines()]
    GSE.close()

    file = open(map, "r")
    GSEtoGSM = defaultdict(list)
    for line in file.readlines():
        line = line.rstrip()
        gse, gsm = line.split(" ")
        GSEtoGSM[gse].append(gsm)
    file.close()

    output = open("related_GSM_IDs.csv", "w")
    GSMs = []
    for series in GSEs:
        GSMs += GSEtoGSM[series]

    GSMs = set(GSMs)
    for id in GSMs:
        output.write(id+"\n")

    output.close()
    print len(GSMs)
getRelatedGSMs()



