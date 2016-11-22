import os
import re
#
# import re
#
# failedGSESoftDownloads = set()
#
#
# cwd = os.getcwd() + '/QuickXMLs/'
#
#
# for filename in os.listdir(cwd):
#     # GSEName = filename[:filename.find(".txt")]
#
#     with open(cwd + filename, "r") as file:
#         n = 0
#         content = file.read()
#
#         if not content.startswith("^SAMPLE"):
#             failedGSESoftDownloads.add(filename[:filename.find(".xml")])
#
# with open("failedGSMSoftDownloads.txt", "w") as file:
#     for gseid in failedGSESoftDownloads:
#         file.write(gseid + "\n")

full = set()
missing = set()

with open("ALLGSMsID.txt", "r") as file:
    for line in file.readlines():
        line = line.strip()
        full.add(line)

with open("failedGSMSoftDownloads.txt", "r") as file:
    for line in file.readlines():
        line = line.strip()
        missing.add(line)

need = set()
with open("ALLGSMsIDSOFT.txt", "r") as file:
    for line in file.readlines():
        line = line.strip()
        if line not in full:
            need.add(line)

print len(need)

with open("needGSMs.txt", "w") as file:
    for s in need:
        file.write(s+"\n")
    for a in missing:
        file.write(a + "\n")
