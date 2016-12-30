

class GSM:
    def __init__(self, GSMID):
        self.id = GSMID
        self.series = []
        self.platForm = None
        self.characteristics = {}
        self.supplementaryData = {}
        self.relations = {}
        self.libraryStrategy = None
        self.SRA = None
        self.type = None
        self.features = None
        self.title = None
        self.InstrumentID = None
        self.organism = None

    # potential attributes
        self.antibody = {}
        self.treatment = {}
        self.tissue = None
        self.disease = None
        self.cellLine = None
        self.cellType = None
        self.genotype = {}
        self.title_found = False
        self.ab_found = False
        self.title_ab = False
        self.input = None
