from ftplib import FTP
import urllib
import os
import tarfile
import re
import argparse
from collections import defaultdict
import csv

def getConnection(ftpAddress='ftp.ncbi.nlm.nih.gov', user='anonymous', passwd='bxia@houstonmethodist.org'):
    ftp = FTP(ftpAddress)
    ftp.login(user, passwd)
    return ftp

def DownloadXml(ftp, GEO_prefix='/geo/series/'):
    files_FTP_surfix = []
    ftp.cwd(GEO_prefix)
    GSE_nnnList = ftp.nlst()
    for GSE_nnn in GSE_nnnList:
        cwd = GEO_prefix+GSE_nnn
        ftp.cwd(cwd)
        GSE_SeriesID = ftp.nlst()
        files_FTP_surfix = files_FTP_surfix + GSE_SeriesID
    return files_FTP_surfix


def searchForData(seriesIDs, feature='h3k4me3',
                  prefix='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
                  outputPrefix="GSEIDwith",
                  ignorecase=re.IGNORECASE):
    IDWithFeatures = []
    GSMs = defaultdict(set)
    failed = []
    for seriesID in seriesIDs:
        # try:
        seriesID = seriesID.rstrip()
        url = prefix + seriesID.rstrip()

        # url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM175628'

        page = urllib.urlopen(url).read()

        # ids = re.search(feature, page, flags=ignorecase)

        gsmids = re.findall("GSM[0-9]+", page)

        for id in gsmids:
            GSMs[id].add(seriesID)



        # except:
        #     failed.append(seriesID)

    with open("additionalGSEtoGSM", "w") as file:
        writer = csv.writer(file)
        for key, value in GSMs.items():
            writer.writerow([key, value])
    # with open("failed"+feature+seriesIDs[0], "w") as file:
    #     for item in failed:
    #         file.write("%s\n" % item)
    return GSMs

def closeConnection(ftp):
    print "job done!"
    ftp.quit()

def py_files(members):
    for tarinfo in members:
        if os.path.splitext(tarinfo.name)[-1] == ".xml":
            yield tarinfo

def openFileGetList(directories):
    results = []
    if isinstance(directories, str):
        directories = [directories]
    for directory in directories:
        with open(directory, "r") as file:
            for line in file.readlines():
                results.append(line.rstrip())
    return results

def getFileList(path=None, partName=None):
    if path == None:
        return []
    elif partName == None:
        return os.listdir(path)
    results = []
    listFiles = os.listdir(path)
    for name in listFiles:
        if name.find(partName) != -1:
            results.append(name)
    return results


#
# files = getFileList(os.getcwd(), 'ALLGSMGSM')
# list = openFileGetList(files)
list = ['GSE59', 'GSE37', 'GSE36', 'GSE35', 'GSE34', 'GSE33', 'GSE32', 'GSE31', 'GSE30', 'GSE39', 'GSE38', 'GSE2', 'GSE5', 'GSE4', 'GSE23149', 'GSE6', 'GSE9', 'GSE8', 'GSE36030', 'GSE70082', 'GSE70083', 'GSE70081', 'GSE62483', 'GSE20', 'GSE21', 'GSE22', 'GSE23', 'GSE24', 'GSE25', 'GSE26', 'GSE27', 'GSE28', 'GSE29', 'GSE59461', 'GSE78897', 'GSE79823', 'GSE79821', 'GSE27464', 'GSE45457', 'GSE83252', 'GSE34554', 'GSE85156', 'GSE57828', 'GSE11062', 'GSE62', 'GSE66866', 'GSE55', 'GSE54', 'GSE57', 'GSE56', 'GSE51', 'GSE50', 'GSE53', 'GSE52', 'GSE66530', 'GSE85', 'GSE64195', 'GSE58', 'GSE18264', 'GSE18262', 'GSE86768', 'GSE19561', 'GSE51564', 'GSE47', 'GSE44', 'GSE45', 'GSE42', 'GSE43', 'GSE40', 'GSE41', 'GSE48', 'GSE68457', 'GSE86718', 'GSE78978', 'GSE79', 'GSE78', 'GSE86675', 'GSE86678', 'GSE70440', 'GSE22441', 'GSE22440', 'GSE43915', 'GSE34061', 'GSE85467', 'GSE73', 'GSE72', 'GSE71', 'GSE70', 'GSE77', 'GSE76', 'GSE75', 'GSE74', 'GSE80788', 'GSE80789', 'GSE11077', 'GSE57814', 'GSE24779', 'GSE46', 'GSE7', 'GSE70099', 'GSE71854', 'GSE24777', 'GSE22438', 'GSE85157', 'GSE68', 'GSE69', 'GSE64', 'GSE65', 'GSE66', 'GSE67', 'GSE60', 'GSE61', 'GSE75898', 'GSE63', 'GSE44288', 'GSE44286', 'GSE59591', 'GSE19654', 'GSE60626', 'GSE35268', 'GSE31477', 'GSE49', 'GSE35774', 'GSE19', 'GSE18', 'GSE11', 'GSE76153', 'GSE13', 'GSE12', 'GSE15', 'GSE14', 'GSE17', 'GSE16', 'GSE40918', 'GSE78129', 'GSE86706', 'GSE82', 'GSE83', 'GSE80', 'GSE81', 'GSE86', 'GSE87', 'GSE84', 'GSE60211', 'GSE89', 'GSE10', 'GSE34691', 'GSE83249', 'GSE33347', 'GSE20512', 'GSE86715', 'GSE50488', 'GSE86712']

GSMlist = searchForData(list)
#

print GSMlist

# GSMwithFeature = searchForData(GSMlist)[0]
# GSEwithFeature = searchForData(GSMlist)[0]
#
# print GSMwithFeature, len(GSMwithFeature)
# print GSEwithFeature, len(GSEwithFeature)



#
# ftp = getConnection()
# list = DownloadXml(ftp)
# # print list
# with open('C:/Users/grunt/Desktop/XMLGEO/Downloads.txt', 'w') as file:
#     for item in list:
#         file.write("%s\n" % item)



# parser = argparse.ArgumentParser()
# parser.add_argument('filename', help='a file contains GSE')
#
# args = parser.parse_args()
#
#
# searchForData(args.filename)




# print grabFiles(ftp, list)
# closeConnection(ftp)





