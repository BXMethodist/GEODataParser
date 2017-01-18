

f1 = open("./uniqueGSM_GEOsearch.txt", "r")

set1= set()

for line in f1.readlines():
    set1.add(line.strip())

f2 = open("./unique_result.txt", "r")

set2 = set()

for line in f2.readlines():
    set2.add(line.strip())


for i in set1:
    if i not in set2:
        print i