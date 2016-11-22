import csv
import re

def GEOSearchResultParser(path=None):
    if path == None:
        return "No path input!"

    uniqueGSM = set()
    with open(path, 'r') as file:
        for line in file.readlines():
            GSEs = re.findall('Accession: GSM[0-9]+', line, flags=0)
            for GSE in GSEs:
                uniqueGSM.add(GSE[11:])

    with open("uniqueGSM.txt", 'w') as file:
        for GSE in uniqueGSM:
            file.write(GSE+'\n')
    return uniqueGSM


list = GEOSearchResultParser("gds_result.txt")

print len(list)
print list

