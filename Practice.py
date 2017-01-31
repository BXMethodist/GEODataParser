import os, re


file = open("gds_result_title.txt" , "r")

info = file.readlines()

result = []

for line in info:
    result +=re.findall('\sGSM[0-9]+\s', line)

output = open("title_gsm.txt", "w")

for r in result:
    output.write(r+"\n")

