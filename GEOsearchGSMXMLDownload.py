from ftplib import FTP
import urllib
import os
import tarfile
import argparse


def grabFiles(directories, urlprefix='/geo/series/',
              file_surfix='.txt', partitionNumber=0):
    failDownloads = []
    for directory in directories:
        downloadAddress = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+directory+"&targ=self&form=text&view=quick"
        print "downloading " + directory + " from " + downloadAddress
        # try:
        urllib.urlretrieve(downloadAddress, "/Users/boxia/Desktop/PycharmProjects/GEODataParser/GSESoftQuick/"+directory+file_surfix)
        # except:
        #
        #     error = 1
        #     while error <= 5:
        #         try:
        #             urllib.urlretrieve(downloadAddress, './Users/boxia/Desktop/PycharmProjects/GEODataParser/GSESoftQuick/'+directory+file_surfix)
        #             break
        #         except:
        #             error += 1
        #     if error >= 5:
        #         failDownloads.append(downloadAddress)
        #         print error
        #     continue
    print failDownloads
    with open('failDownLoadGSE'+str(partitionNumber), 'w') as file:
        for item in failDownloads:
            file.write("%s\n" % item)

def getPartition(length, number):
    numberSamplesPartition = length / number
    numberSamplesRemainder = length % number
    partitions = [
        [start * numberSamplesPartition, (start + 1) * numberSamplesPartition] if start != number - 1 else [
            start * numberSamplesPartition, (start + 1) * numberSamplesPartition + numberSamplesRemainder] for start in
        range(number)]
    return partitions

def MinimlDowloader(numberPartition, numberofPartitions=1):
    numberPartition = int(numberPartition)
    numberofPartitions = int(numberofPartitions)
    GSMs = []
    with open("/Users/boxia/Desktop/PycharmProjects/GEODataParser/failedGSESoftDownloads.txt", "r") as file:
        for line in file.readlines():
            GSMs.append(line.rstrip())
    start, end = getPartition(len(GSMs), numberofPartitions)[numberPartition]

    downloadFiles = GSMs[start: end]
    grabFiles(downloadFiles, partitionNumber=numberPartition)

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('number', help='a file contains GSE', type=int)
#
#     args = parser.parse_args()
#
#     MinimlDowloader(args.number)


MinimlDowloader(0)


