from GSM import *
import re
import os

def SoftFileParser(cwd=None):
    if cwd == None:
        cwd = os.getcwd()+'\SOFT'

    samples = {}
    for filename in os.listdir(cwd):
        with open(filename, "r") as f:
            line = f.readline()
            while True:
                if not line:
                    break
                if line.startswith("^SERIES"):
                    info = line.rstrip('\n').split(" ")
                    series_number = info[-1]
                if line.startswith("!Series_sample_id"):
                    info = line.rstrip('\n').split(" ")
                    sample_id = info[-1]
                    samples[sample_id] = GSM(sample_id, series_number)
                if line.startswith("^SAMPLE"):
                    info = line.rstrip('\n').split(" ")
                    sample_id = info[-1]
                    sample = samples[sample_id]
                    while True:
                        line = f.readline()
                        if not line:
                            break
                        if line.startswith("^SAMPLE"):
                            break
                        if line.startswith("!Sample_characteristics_ch1"):
                            info = re.split("=|:", line.rstrip('\n'))
                            sample.characteristics[info[1].strip().upper()] = info[-1]
                    continue
                line = f.readline()
    return samples
