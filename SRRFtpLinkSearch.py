from ftplib import FTP


def getConnection(ftpAddress='ftp-trace.ncbi.nlm.nih.gov', user='anonymous', passwd='bxia@houstonmethodist.org'):
    ftp = FTP(ftpAddress)
    ftp.login(user, passwd)
    return ftp


def getDownloadLinks(ftp, HttplinksFile):
    ftplinks = []

    Httplinks = []

    failed = []

    with open(HttplinksFile, "r") as file:
        for line in file.readlines():
            Httplinks.append(line.strip())

    Httplinks = ["https://www.ncbi.nlm.nih.gov/sra?term=SRX532969"]

    for link in Httplinks:
        # try:
        SRX = link[link.find("SRX"):].strip()
        SRXfolder = SRX[:6]

        subfolderAddress = "/sra/sra-instant/reads/ByExp/sra/SRX/"+SRXfolder+'/'+SRX+'/'

        ftp.cwd(subfolderAddress)

        SRRfolderList = ftp.nlst()

        for SRRfolder in SRRfolderList:
            SRRfolderAddress = subfolderAddress + SRRfolder + '/'
            ftp.cwd(SRRfolderAddress)

            SRRfiles = ftp.nlst()
            for SRRfile in SRRfiles:
                downloadFTP = ftp.host  + SRRfolderAddress + SRRfile
                ftplinks.append(downloadFTP)
        # except:
        #     failed.append(SRX)

    with open("HumanWithH3K4me3Links.txt", "w") as file:
        for link in ftplinks:
            file.write("%s\n" %link)

    with open("failedSRX.txt", "w") as file:
        for srx in failed:
            file.write("%s\n" % srx)

    return ftplinks

ftp = getConnection()

getDownloadLinks(ftp, "HumanWithH3K4me3Download.txt")