import urllib2, re
from collections import defaultdict
import psutil

def search_term_to_GSM(terms):
    result_ids = set()
    result_gsms = set()

    for term in terms:
        cur_result_ids = readid("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GSM%5BETYP%5D"+term+"&retmax=10000")

        result_ids = result_ids.union(cur_result_ids)

    result_ids = list(result_ids)

    for i in range(0, len(result_ids), 500):
        ids = ",".join(result_ids[i:i+500])
        gsm_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gds&version=2.0&id=" + ids

        f = urllib2.urlopen(gsm_url)
        content = f.read()
        result_gsms = result_gsms.union(re.findall('\sGSM[0-9]+\s', content))
    gsms = set()

    for gsm in result_gsms:
        gsms.add(gsm.strip())

    return gsms


def readid(url):
    f = urllib2.urlopen(url)

    content = f.read()
    content = content[content.find("<IdList>") + 8:content.find("</IdList>")-1]
    content = content.replace("<Id>", " ")
    content = content.replace("</Id>", " ")
    result_ids = set(content.split())
    return result_ids


class GSM:
    def __init__(self, GSM_id):
        self.id = GSM_id
        self.series = []
        self.platForm = ""
        self.characteristics = {}
        self.supplementaryData = {}
        self.relations = {}
        self.libraryStrategy = ""
        self.SRA = ""
        self.type = ""
        self.features = ""
        self.title = ""
        self.InstrumentID = ""
        self.organism = ""
        self.url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="+GSM_id+"&targ=self&form=text&view=quick"

    # potential attributes
        self.antibody = {}
        self.treatment = {}
        self.tissue = ""
        self.disease = ""
        self.cellLine = ""
        self.cellType = ""
        self.genotype = {}
        self.title_found = False
        self.ab_found = False
        self.title_ab = False
        self.input = ""
        self.encode = ""
        self.roadmap = ""