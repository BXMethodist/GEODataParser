from ftplib import FTP
import urllib
import os

def getConnection(ftpAddress='ftp.ncbi.nlm.nih.gov', user='anonymous', passwd='bxia@houstonmethodist.org'):
    ftp = FTP(ftpAddress)
    ftp.login(user, passwd)
    return ftp

def DownloadSoft(ftp, GEO_prefix='/geo/series/'):
    files_FTP_surfix = []
    ftp.cwd(GEO_prefix)
    GSE_nnnList = ftp.nlst()
    for GSE_nnn in GSE_nnnList:
        cwd = GEO_prefix+GSE_nnn
        ftp.cwd(cwd)
        GSE_xxxxxList = ftp.nlst()
        GSE_SeriesID = [cwd+'/'+x for x in GSE_xxxxxList]
        print GSE_SeriesID
        files_FTP_surfix = files_FTP_surfix + GSE_SeriesID
    return files_FTP_surfix

def grabFiles(ftp, directories, directory_surfix='/soft/', file_surfix='_family.soft.gz'):
    ftpAddress = ftp.host
    try:
        for directory in directories:
            lastFowardSlash = directory.rfind("/")
            filename = directory[lastFowardSlash + 1:] + file_surfix
            downloadAddress = 'ftp://'+ftpAddress + directory + directory_surfix + filename
            print "downloading " + filename + " from " + downloadAddress
            urllib.urlretrieve(downloadAddress, 'C:/Users/grunt/Desktop/SoftGEO/'+filename)
        return 0
    except:
        raise

def closeConnection(ftp):
    print "job done!"
    ftp.quit()


ftp = getConnection()
list = DownloadSoft(ftp)
grabFiles(ftp, list)
closeConnection(ftp)

