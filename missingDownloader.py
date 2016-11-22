from ftplib import FTP
import urllib
import os
import tarfile
import argparse

def getConnection(ftpAddress='ftp.ncbi.nlm.nih.gov', user='anonymous', passwd='bxia@houstonmethodist.org'):
    ftp = FTP(ftpAddress)
    ftp.login(user, passwd)
    return ftp


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