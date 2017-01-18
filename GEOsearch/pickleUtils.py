import pickle, os
import sqlite3
from collections import defaultdict
from ftplib import FTP
import pandas as pd

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


def getConnection(ftpAddress='ftp-trace.ncbi.nlm.nih.gov', user='anonymous', passwd='bxia@houstonmethodist.org'):
    ftp = FTP(ftpAddress)
    ftp.login(user, passwd)
    return ftp


def GSMGSESRR_pickle(path="/home/tmhbxx3/scratch/XMLhttp/GSESoftQuick/"):
    GSMGSE_map = defaultdict(set)

    import psutil
    proc = psutil.Process()

    gse_files = os.listdir(path)

    for gse_file in gse_files:
        gse_file_obj = open(path + gse_file, "r")
        info = gse_file_obj.readlines()
        gse_file_obj.close()

        gse_id = gse_file[:-4]

        if len(proc.open_files()) > 100:
            print "More than one file is open, stop!"
            print proc.open_files()
            return

        for line in info:
            if line.startswith("!Series_sample_id"):
                gsm = line.split("=")[1].strip()
                GSMGSE_map[gsm].add(gse_id)
                GSMGSE_map[gse_id].add(gsm)

    save_obj(GSMGSE_map, "GSMGSE_map")
    return


def roadmap_encode(path="/home/tmhbxx3/scratch/XMLhttp/GSESoftQuick/"):
    roadmap = set()
    encode = set()
    import psutil
    proc = psutil.Process()

    gse_files = os.listdir(path)

    for gse_file in gse_files:
        gse_file_obj = open(path + gse_file, "r")
        info = gse_file_obj.readlines()
        gse_file_obj.close()

        if len(proc.open_files()) > 100:
            print "More than one file is open, stop!"
            print proc.open_files()
            return

        for line in info:
            if line.startswith("!Series_project"):
                if line.find("Roadmap Epigenomics") != -1:
                    roadmap.add(gse_file[:-4])
                if line.find("ENCODE") != -1:
                    encode.add(gse_file[:-4])

    save_obj(roadmap, "Roadmap_gse")
    save_obj(encode, "ENCODE_gse")
    return
