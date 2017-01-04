import pickle
import sqlite3
from collections import defaultdict
import urllib
import re
from ftplib import FTP

def save_obj(obj, name):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def getConnection(ftpAddress='ftp-trace.ncbi.nlm.nih.gov', user='anonymous', passwd='bxia@houstonmethodist.org'):
    ftp = FTP(ftpAddress)
    ftp.login(user, passwd)
    return ftp


def GSMGSE_pickle(path='/home/tmhbxx3/archive/GEO_MetaDatabase/geoMetaData.db', partition=None):
    ftp = getConnection()

    db = sqlite3.connect(path)
    db.text_factory = str

    # query = db.execute("select distinct GSM_ID, GSE_ID from GSEtoGSM").fetchall()
    #
    # GSMGSE_map = defaultdict(set)
    # for GSM_ID, GSE_ID in query:
    #     GSMGSE_map[GSM_ID].add(GSE_ID)
    #     GSMGSE_map[GSE_ID].add(GSM_ID)
    #
    # save_obj(GSMGSE_map, "GSMGSE_map")

    query = db.execute("select distinct GSM_ID, SRA from GSM").fetchall()

    GSMSRR_map = defaultdict(set)

    if partition is not None:
        length = len(query)
        blocksize = length/100
        start = partition * blocksize
        end = (partition+1) * blocksize
        query = query[start:end]

    for GSM_ID, SRXlink in query:
        if SRXlink == "":
            continue

        page = urllib.urlopen(SRXlink).read()

        SRRids = re.findall("SRR[0-9]+", page)
        DRRids = re.findall("DRR[0-9]+", page)
        ERXids = re.findall("ERX[0-9]+", page)

        SRRids += DRRids + ERXids

        GSMSRR_map[GSM_ID].union(set(SRRids))

        for ID in SRRids:
            GSMSRR_map[ID].add(GSM_ID)

    if partition is not None:
        output_name_GSMSRR = "GSMSRR_map"+str(partition)
        output_name_SRR_ftp_map = "SRR_ftp" + str(partition)
    else:
        output_name_GSMSRR = "GSMSRR_map"
        output_name_SRR_ftp_map = "SRR_ftp"
    save_obj(GSMSRR_map, output_name_GSMSRR)


    # SRR_ftp_map = defaultdict(set)
    # save_obj(SRR_ftp_map, output_name_SRR_ftp_map)

    db.close

if __name__ == "__main__":
    GSMGSE_pickle(partition=None)