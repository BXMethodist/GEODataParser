"""
Copyright (c) <2017> <Dr. Kaifu Chen lab, Research Institute, Houston Methodist Hospital >

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from ftplib import FTP
from setup import get_settings
from pickleUtils import load_obj, save_obj
import pandas as pd, urllib2, requests, mimetypes
import gc, sqlite3, json


def connectToGEO(user='anonymous', ftpAddress='ftp.ncbi.nlm.nih.gov'):
    ### create GEO ftp connection to NCBI
    parameters = get_settings()
    email = parameters['email']

    ftp = FTP(ftpAddress)
    ftp.login(user, email)
    return ftp


def updateGSMGSE_Encode_Roadmap(GSMGSE_map, Encode_map, Roadmap_map, GGR_map, MetaData_path):
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

    if MetaData_path is not None:
        db = sqlite3.connect(MetaData_path)
    else:
        db = None

    for id in seriesID:
        if id not in local_GSEs:
            result = GSE_info(id)
            if result is None:
                continue
            newGSMs, encode, roadmap, ggr = result

            newGSMs_info = []
            curGSMGSE_map = {}
            for gsm in newGSMs:
                related_GSEs = GSM_info(gsm)
                if related_GSEs is not None:
                    curGSMGSE_map[gsm] = related_GSEs
                    if db is not None:
                        info = downloadGSM(gsm)
                        if info is not None:
                            try:
                                metadata = json.dumps(info)
                            except:
                                new_info = []
                                for i in info:
                                    new_info.append(unicode(i, errors='ignore'))
                                metadata = json.dumps(new_info)
                            newGSMs_info.append((gsm, metadata))

            if len(newGSMs_info) != len(newGSMs):
                print id, ", try to update, but failed to get all gsms"
            else:
                GSMGSE_map[id] = newGSMs
                if encode:
                    Encode_map.add(id)
                if roadmap:
                    Roadmap_map.add(id)
                if ggr:
                    GGR_map.add(id)
                for gsm in newGSMs:
                    GSM_need_update.add(gsm)
                    GSMGSE_map[gsm] = curGSMGSE_map[gsm]
                for data in newGSMs_info:
                    gsm, metadata = data
                    db.execute("insert into GSM values(?, ?)", (gsm, metadata))
                    print "update ", gsm
    db.commit()
    db.close()
    return GSMGSE_map, Encode_map, Roadmap_map, GGR_map, GSM_need_update


def updateGSMSRR(GSMSRR_map, GSMs):
    for gsm in GSMs:
        try:
            df = pd.read_csv(
                "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + gsm,
                index_col=0)

            SRRs = set(df.index.values)
            GSMSRR_map[gsm] = SRRs

            for srr in SRRs:
                GSMSRR_map[srr] = gsm
        except:
            pass

    return GSMSRR_map


def GSE_info(id):
    try:
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + id + "&targ=self&form=text&view=quick"
        encode = False
        roadmap = False
        ggr = False
        gsms = set()

        web = urllib2.urlopen(url)
        info = web.readlines()
        web.close()

        for line in info:
            if line.startswith("!Series_project"):
                if line.find("Roadmap Epigenomics") != -1:
                    roadmap = True
                if line.find("ENCODE") != -1:
                    encode = True
                if line.find("GGR") != -1:
                    ggr = True
            if line.startswith("!Series_sample_id"):
                gsm = line.split("=")[1].strip()
                gsms.add(gsm)
        return gsms, encode, roadmap, ggr
    except:
        print 'failed to get metadata for GSE', ' ', id
        return None


def GSM_info(id):
    try:
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + id + "&targ=self&form=text&view=quick"

        web = urllib2.urlopen(url)
        info = web.readlines()
        web.close()

        gses = set()
        for line in info:
            if line.startswith("!Sample_series_id"):
                gses.add(line[line.find("=")+1:].strip())
        return gses
    except:
        print 'failed to get GSE id for ', id
        return None


def downloadGSM(gsmID):
    ## download GSM meta data from NCBI
    try:
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + gsmID + "&targ=self&form=text&view=quick"
        response = requests.get(url)
        content_type = response.headers['content-type']
        extension = mimetypes.guess_extension(content_type)
        if content_type == "geo/text" and extension != ".html":
            web = urllib2.urlopen(url)
            info = web.readlines()
            web.close()
            del web
            gc.collect()
        else:
            info = None
        return info
    except:
        print 'failed to download metadata of ', gsmID
        return None

def update():
    parameters = get_settings()
    GSMGSE_map = load_obj(parameters['GSMGSE_pkl_path'])
    GSMSRR_map = load_obj(parameters['GSMtoSRRpkl'])
    Encode_map = load_obj(parameters['Encode'])
    Roadmap_map = load_obj(parameters['Roadmap'])
    GGR_map = load_obj(parameters['GGR'])

    MetaData_path = parameters["MetaData"]

    if MetaData_path == "None":
        MetaData_path = None

    GSMGSE_map, Encode_map, Roadmap_map, GGR_map, GSM_need_update = \
        updateGSMGSE_Encode_Roadmap(GSMGSE_map, Encode_map, Roadmap_map, GGR_map, MetaData_path)

    GSMSRR_map = updateGSMSRR(GSMSRR_map, GSM_need_update)

    save_obj(GSMGSE_map, parameters['GSMGSE_pkl_path'][:-4])
    save_obj(Encode_map, parameters['Encode'][:-4])
    save_obj(Roadmap_map, parameters['Roadmap'][:-4])
    save_obj(GGR_map, parameters['GGR'][:-4])
    save_obj(GSMSRR_map, parameters['GSMtoSRRpkl'][:-4])




