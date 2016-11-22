from ftplib import FTP
import urllib
import os
import tarfile
import argparse

def getConnection(ftpAddress='ftp.ncbi.nlm.nih.gov', user='anonymous', passwd='bxia@houstonmethodist.org'):
    ftp = FTP(ftpAddress)
    ftp.login(user, passwd)
    return ftp

def DownloadXml(ftp, outputPath, GEO_prefix='/geo/series/', numberPartition=1):
    files_FTP_surfix = []
    seriesID = []
    print ftp.host
    ftp.cwd(GEO_prefix)
    GSE_nnnList = ftp.nlst()
    for GSE_nnn in GSE_nnnList:
        cwd = GEO_prefix+GSE_nnn
        ftp.cwd(cwd)
        GSE_xxxxxList = ftp.nlst()
        seriesID += [id for id in GSE_xxxxxList]
        GSE_SeriesID = [cwd+'/'+x for x in GSE_xxxxxList]

        files_FTP_surfix = files_FTP_surfix + GSE_SeriesID

    with open(outputPath, 'w') as file:
        for item in files_FTP_surfix:
            file.write("%s\n" % item)

    numberSamplesPartition = len(files_FTP_surfix)/numberPartition
    numberSamplesRemainder = len(files_FTP_surfix)%numberPartition
    partitions = [[start*numberSamplesPartition, (start+1)*numberSamplesPartition] if start != numberPartition-1 else [start*numberSamplesPartition, (start+1)*numberSamplesPartition+numberSamplesRemainder] for start in range(numberPartition)]

    with open("GSEIDs", 'w') as file:
        for item in seriesID:
            input = item+ "\n"
            file.write(input)

    return files_FTP_surfix, partitions

def openFileGetList(directories):
    GSEs = []
    with open(directories, "r") as file:
        for line in file.readlines():
            GSEs.append(line.rstrip())
    return GSEs

def grabFiles(ftp, directories, urlprefix='/geo/series/',
              directory_surfix='/miniml/', file_surfix='_family.xml.tgz', partitionNumber=0):
    ftpAddress = ftp.host
    failDownloads = []
    for directory in directories:
        index = directory.rfind('/')+1
        GSEID= directory[index:]
        filename = GSEID+file_surfix
        downloadAddress = 'ftp://'+ftpAddress+directory+directory_surfix + filename
        print "downloading " + filename + " from " + downloadAddress
        try:
            urllib.urlretrieve(downloadAddress, './MINIMLs/'+filename)
        except:

            error = 1
            while error <= 5:
                try:
                    urllib.urlretrieve(downloadAddress, '/MINIMLs/' + filename)
                    break
                except:
                    error += 1
            failDownloads.append(downloadAddress)
            continue
        tar = tarfile.open('./MINIMLs/'+filename)
        tar.extractall(path='./MINIMLs/',members=py_files(tar))
        tar.close()
        os.remove('./MINIMLs/'+filename)
    with open('failDownloadsMINIMLs'+str(partitionNumber), 'w') as file:
        for item in failDownloads:
            file.write("%s\n" % item)

def closeConnection(ftp):
    print "job done!"
    ftp.quit()

def py_files(members):
    for tarinfo in members:
        if os.path.splitext(tarinfo.name)[-1] == ".xml":
            yield tarinfo


def MinimlDowloader(ftpGSElist, outputPath=".", numberPartition=0, numberofPartitions=1):
    ftp = getConnection()
    numberPartition = int(numberPartition)
    numberofPartitions = int(numberofPartitions)
    GSEsPrefix, partitionFile = DownloadXml(ftp, ftpGSElist, numberPartition=numberofPartitions)
    print len(GSEsPrefix)
    # ftpPrefixes = openFileGetList(ftpPrefixList)
    #
    # partitions = []
    # with open(ftpPrefixList+"partition", "r") as file:
    #     for line in file.readlines():
    #         partitions.append(line.rstrip().split(" "))


    # downloadFiles = GSEsPrefix[partitionFile[numberPartition][0]: partitionFile[numberPartition][1]]
    # grabFiles(ftp, downloadFiles, partitionNumber=numberPartition)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('filename', help='a file contains GSE')
#     # parser.add_argument('numberofPartitions', help='how many partitions')
#     # parser.add_argument('numberofPartition', help='which partition')
#
#     args = parser.parse_args()
#
#     MinimlDowloader(args.filename)





#
MinimlDowloader("ftplist1", 0, 1)
