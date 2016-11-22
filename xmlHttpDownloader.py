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

    # with open(outputPath, 'w') as file:
    #     for item in files_FTP_surfix:
    #         file.write("%s\n" % item)

    numberSamplesPartition = len(files_FTP_surfix)/numberPartition
    numberSamplesRemainder = len(files_FTP_surfix)%numberPartition
    partitions = [[start*numberSamplesPartition, (start+1)*numberSamplesPartition] if start != numberPartition-1 else [start*numberSamplesPartition, (start+1)*numberSamplesPartition+numberSamplesRemainder] for start in range(numberPartition)]

    # with open(outputPath+'partition', 'w') as file:
    #     for item in partitions:
    #         input = str(item[0])+" "+str(item[1])+ "\n"
    #         file.write(input)

    return files_FTP_surfix, partitions, seriesID

def grabFiles(directories, urlprefix='/geo/series/',
              file_surfix='.xml', partitionNumber=0):
    failDownloads = []
    for directory in directories:
        downloadAddress = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+directory+"&targ=self&form=xml&view=quick"
        # print "downloading " + directory + " from " + downloadAddress
        # try:
        urllib.urlretrieve(downloadAddress, './XMLs/'+directory+file_surfix)
        # except:

            # error = 1
            # while error <= 5:
            #     try:
            #         urllib.urlretrieve(downloadAddress, './XMLs/'+directory+file_surfix)
            #         break
            #     except:
            #         error += 1
            # failDownloads.append(downloadAddress)
            # continue
    #     tar = tarfile.open('./MINIMLs/'+filename)
    #     tar.extractall(path='./MINIMLs/',members=py_files(tar))
    #     tar.close()
    #     os.remove('./MINIMLs/'+filename)
    # with open('failDownloadsMINIMLs'+str(partitionNumber), 'w') as file:
    #     for item in failDownloads:
    #         file.write("%s\n" % item)

def closeConnection(ftp):
    print "job done!"
    ftp.quit()

def py_files(members):
    for tarinfo in members:
        if os.path.splitext(tarinfo.name)[-1] == ".xml":
            yield tarinfo


def MinimlDowloader(numberPartition, numberofPartitions=10):
    ftpGSElist = "GSEID"
    ftp = getConnection()
    numberPartition = int(numberPartition)
    numberofPartitions = int(numberofPartitions)
    GSEs = DownloadXml(ftp, ftpGSElist, numberPartition=numberofPartitions)[2]
    partitionFile = DownloadXml(ftp, ftpGSElist, numberPartition=numberofPartitions)[1]

    downloadFiles = GSEs[partitionFile[numberPartition][0]: partitionFile[numberPartition][1]]
    grabFiles(downloadFiles, partitionNumber=numberPartition)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('number', help='a file contains GSE', type=int)
    # parser.add_argument('numberofPartitions', help='how many partitions')
    # parser.add_argument('numberofPartition', help='which partition')

    args = parser.parse_args()

    MinimlDowloader(args.number)





#
# MinimlDowloader("ftplist1", 1, 10)
