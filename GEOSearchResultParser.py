import csv
import re

def GEOSearchResultParser(path=None):
    if path == None:
        return "No path input!"

    uniqueGSM = []
    with open(path, 'r') as file:
        for line in file.readlines():
            GSMs = re.findall('Accession: GSM[0-9]+', line, flags=0)
            for GSM in GSMs:
                uniqueGSM.append(GSM[11:])

    with open("uniqueGSM_GEOsearch.txt", 'w') as file:
        for GSE in uniqueGSM:
            file.write(GSE+'\n')
    return uniqueGSM


list = GEOSearchResultParser("gds_result_12_12_for_database_sync.txt")

print len(list)
print list

