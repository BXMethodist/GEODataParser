from subprocess import call
import os

file = open("./wgetFiles11-22.txt", "r")

listFiles = os.listdir("/home/tmhbxx3/archive/H3K4me3")
for line in file.readlines():
    fileName = line[line.rfind("SRR"):].rstrip()
    if fileName not in listFiles:
        commands = line.split(" ")
        call(commands)

file.close()
