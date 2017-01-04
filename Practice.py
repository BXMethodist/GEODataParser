import psutil

file_obj = open("./csv/humanWithH3K4me3.csv", "r")

file_obj2 = open("./csv/sampleWithH3K4me3.csv", "r")

proc = psutil.Process()
print proc.open_files()
info = file_obj.readlines()
file_obj.close()
file_obj2.close()
print len(info)
print proc.open_files()

# for line in info:
#     print line

