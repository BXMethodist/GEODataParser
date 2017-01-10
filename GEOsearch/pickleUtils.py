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

    print "fetching complete!"

    GSMSRR_map = defaultdict(set)
    SRR_map = {}
    failed = []

    if partition is not None:
        length = len(query)
        blocksize = length/1000
        start = partition * blocksize
        end = (partition+1) * blocksize
        query = query[start:end]

    for GSM_ID, SRXlink in query:
        if SRXlink == "":
            continue
        try:
            SRX = SRXlink[(SRXlink.find("=") + 1):].strip()
            url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?sp=runinfo&acc=" + SRX + "&retmode=txt"
            page = pd.read_csv(url)
            for i in range (page.shape[0]):
                SRRid = page.ix[i, "Run"]
                download_path = page.ix[i, "download_path"]
                layout = page.ix[i, "LibraryLayout"]
                avg_length = page.ix[i, "avgLength"]
                GSMSRR_map[GSM_ID].add(SRRid)
                SRR_map[SRRid] = SRR(SRRid, download_path, layout, avg_length, GSM_ID)
        except:
            failed.append(SRXlink)

    if partition is not None:
        output_name_GSMSRR = "GSMSRR_map"+str(partition)
        output_name_SRR_ftp_map = "SRR_map" + str(partition)
    else:
        output_name_GSMSRR = "GSMSRR_map"
        output_name_SRR_ftp_map = "SRR_map"
    save_obj(GSMSRR_map, output_name_GSMSRR)

    save_obj(SRR_map, output_name_SRR_ftp_map)

    for srx in failed:
        print srx

    db.close


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


class SRR:
    def __init__(self, SRRid, download_path, layout, avg_length, gsm_id):
        self.SRRid = SRRid
        self.download_path = download_path
        self.layout = layout
        self.avg_length = avg_length
        self.gsm_id = gsm_id

# if __name__ == "__main__":
#     # GSMGSE_pickle(partition=None)
#     roadmap_encode()