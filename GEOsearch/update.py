
from ftplib import FTP
from setup import get_settings
from pickleUtils import load_obj, save_obj
import pandas as pd, urllib


def connectToGEO(user='anonymous', ftpAddress='ftp.ncbi.nlm.nih.gov'):
    ### create GEO ftp connection to NCBI
    parameters = get_settings()
    email = parameters['email']

    ftp = FTP(ftpAddress)
    ftp.login(user, email)
    return ftp


def updateGSMGSE_Encode_Roadmap(GSMGSE_map, Encode_map, Roadmap_map):
    # looking for new GSE, if not in local, update the pkl
    local_GSEs = set()
    local_GSMs = set()

    for key in GSMGSE_map.keys():
        if key.startswith("GSM"):
            local_GSMs.add(key)
        else:
            local_GSEs.add(key)

    ftp = connectToGEO()

    GEO_prefix = '/geo/series/'
    seriesID = []
    ftp.cwd(GEO_prefix)
    GSE_nnnList = ftp.nlst()
    for GSE_nnn in GSE_nnnList:
        cwd = GEO_prefix+GSE_nnn
        ftp.cwd(cwd)
        GSE_xxxxxList = ftp.nlst()
        seriesID += [id for id in GSE_xxxxxList]

    seriesID = set(seriesID)

    GSM_need_update = set()

    for id in seriesID:
        if id not in local_GSEs:
            newGSMs, encode, roadmap = GSE_info(id)

            if encode:
                Encode_map.add(id)
            if roadmap:
                Roadmap_map.add(id)
            GSMGSE_map[id] = newGSMs
            for gsm in newGSMs:
                related_GSEs = GSM_info(gsm)
                GSMGSE_map[gsm] = related_GSEs
                GSM_need_update.add(gsm)

    return GSMGSE_map, Encode_map, Roadmap_map, GSM_need_update


def updateGSMSRR(GSMSRR_map, GSMs):
    for gsm in GSMs:
        df = pd.read_csv(
            "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + gsm,
            index_col=0)

        SRRs = set(df.index.values)
        GSMSRR_map[gsm] = SRRs

        for srr in SRRs:
            GSMSRR_map[srr] = gsm

    return GSMSRR_map


def GSE_info(id):
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + id + "&targ=self&form=text&view=quick"
    encode = False
    roadmap = False
    gsms = set()

    for line in urllib.urlopen(url).readlines():
        if line.startswith("!Series_project"):
            if line.find("Roadmap Epigenomics") != -1:
                roadmap = True
            if line.find("ENCODE") != -1:
                encode = True
        if line.startswith("!Series_sample_id"):
            gsm = line.split("=")[1].strip()
            gsms.add(gsm)
    return gsms, encode, roadmap


def GSM_info(id):
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + id + "&targ=self&form=text&view=quick"

    gses = set()
    for line in urllib.urlopen(url).readlines():
        if line.startswith("!Sample_series_id"):
            gses.add(line[line.find("=")+1:].strip())
    return gses

def update():
    parameters = get_settings()
    GSMGSE_map = load_obj(parameters['GSMGSE_pkl_path'])
    GSMSRR_map = load_obj(parameters['GSMtoSRRpkl'])
    Encode_map = load_obj(parameters['Encode'])
    Roadmap_map = load_obj(parameters['Roadmap'])

    GSMGSE_map, Encode_map, Roadmap_map, GSM_need_update = updateGSMGSE_Encode_Roadmap(GSMGSE_map, Encode_map, Roadmap_map)

    GSMSRR_map = updateGSMSRR(GSMSRR_map, GSM_need_update)

    save_obj(GSMGSE_map, parameters['GSMGSE_pkl_path'][:-4])
    save_obj(GSMSRR_map, parameters['GSMtoSRRpkl'][:-4])
    save_obj(Encode_map, parameters['Encode'][:-4])
    save_obj(Roadmap_map, parameters['Roadmap'][:-4])

update()