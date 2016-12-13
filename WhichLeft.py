import os
import re

GSMs = {}
with open("./uniqueGSM_GEOsearch.txt", "r") as file:
    n = 0
    for line in file.readlines():
        line = line.rstrip()
        GSMs[line] = n
        n += 1

cwd = "./GEOSearchXMLs"

files = os.listdir(cwd)
indexes = []
for file in files:
    result = re.findall("GSM[0-9]+", file)
    if result != []:
        indexes.append(GSMs[result[0]])

indexes.sort()
results = []
for i in range(1, len(indexes)):
    if indexes[i]-1 > indexes[i-1] and indexes[i]-2 > indexes[i-1]:
        # print indexes[i]
        results.append(indexes[i-1])
print indexes[-1]
results.append(indexes[-1])
print results