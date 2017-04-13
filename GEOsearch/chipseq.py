import sqlite3, json
from multiprocessing import Process, Queue
from sample import GSM
from input_search_utils import get_MetaInfo
from setup import get_settings
from pickleUtils import load_obj
from collections import defaultdict
import re


def chip_filter(gsms, queue, excludeGSEs):
    db = sqlite3.connect('/home/tmhbxx3/scratch/GCF_OL/pkl/geoMetaData.db')

    results = []

    GSE_GSM = defaultdict(set)

    count = 0

    for sampleName in gsms:
        sample = GSM(sampleName)

        sampleTitle = ""
        sampleType = ""
        sampleLibraryStrategy = ''
        sampleOrganism = ''
        samplePlatForm = ''
        sampleInstrumentID = ''

        characteristics = defaultdict(str)
        supplementaryData = defaultdict(str)
        relations = defaultdict(str)
        sampleSeriesID = set()
        target_feature = {}

        antibody = {}
        treatment = {}
        tissue = ""
        disease = ""
        cellLine = ""
        cellType = ""
        genoType = {}

        info = get_MetaInfo(db, sample, count)

        count += 1

        if info is None:
            continue

        for line in info:
            if line.startswith("!Sample_title"):
                sampleTitle = line[line.find("=") + 1:].strip()
                if sampleTitle.find(";") != -1:
                    sampleTitle = sampleTitle[:sampleTitle.find(";")]
            if line.startswith("!Sample_type"):
                sampleType = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_organism"):
                sampleOrganism = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_characteristics_ch"):
                characteristic = line[line.find("=") + 1:].strip()
                key, value = characteristic[:characteristic.find(":")].strip(), characteristic[
                                                                                characteristic.find(":") + 1:].strip()
                if key in characteristics:
                    characteristics[key] += ", " + value
                else:
                    characteristics[key] = value
            if line.startswith("!Sample_platform_id "):
                samplePlatForm = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_library_strategy"):
                sampleLibraryStrategy = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_supplementary"):
                dataUrl = line[line.find("=") + 1:].strip()
                keyName = "Data" + str(len(supplementaryData) + 1)
                supplementaryData[keyName] = dataUrl
            if line.startswith("!Sample_relation"):
                relation = line[line.find("=") + 1:].strip()
                key, value = relation[:relation.find(":")].strip(), relation[relation.find(":") + 1:].strip()
                relations[key] = value
            if line.startswith("!Sample_series_id"):
                sampleSeriesID.add(line[line.find("=") + 1:].strip())
            if line.startswith("!Sample_instrument_model"):
                sampleInstrumentID = line[line.find("=") + 1:].strip()
            if line.startswith("!Sample_contact_institute = ENCODE DCC"):
                sample.encode = True

        sample.characteristics = characteristics
        sample.supplementaryData = supplementaryData
        sample.title = sampleTitle
        sample.type = sampleType
        sample.libraryStrategy = sampleLibraryStrategy
        sample.organism = sampleOrganism
        sample.SRA = relations["SRA"]
        sample.series = list(sampleSeriesID)

        sample.platForm = samplePlatForm
        sample.features = target_feature
        sample.InstrumentID = sampleInstrumentID

        for key, value in characteristics.items():
            if key.lower() in ["chip antibody", "chip", "antigen", "antibody", "antibodies", "chip antibodies",
                               "antibody name", "antibody target", "target", "antibody/capture",
                               "antibody/vendor/catalog#",
                               "chip ab", "chip antibody1", "chip antibody2", "chip-antibody", "chip_antibody",
                               "chip-seq antibody", "chip_antibodies", "chip-antibodies", "histone mark",
                               "epigenetic feature",
                               "histone modification", "antibody antibodydescription", "chip antibody (epitope/name)",
                               "factor", "chip antibody/mbd affinity column", "chip/dip antibody", "antibody epiptope",
                               "antibody source", 'modification', "antibody (vendor': ' catalog#, or reference)",
                               "experiment", "purification antibody", "antibody/details", "antibody epiptope",
                               "antibody information", "chip antibody / digestive enzyme", "chip antiboy",
                               "ip antibody", "chip antibody target", "modification", "histone", "enrichment procedure",
                               "antibody (vendor': ' catalog#, or reference)",
                               "developmental stage/condition/extract protocol",
                               "antibody source"] \
                    or re.search('antibody epitope|using[\w\s]+antibod|immunoprecipitat', key, flags=re.IGNORECASE):
                if key in antibody:
                    antibody[key] += ", " + value
                else:
                    antibody[key] = value
            if key.lower() in ["treatment", "condition", "activation stimuli", "cell condition", "cell treatment",
                               "cell-treatment", "drug treatment", "stress", "overexpression", "treatment drug",
                               "treatment group"] \
                    or re.search("(?:dsrna|infect|rnai|shrna|sirna|transduc|transfec|agent[#]*[0-9]*|activat)", key,
                                 flags=re.IGNORECASE):
                treatment[key] = value

            if key.lower() in ["tissue", "body part", "body site", "tissue type"]:
                tissue = value

            if key.lower() in ["cancer type", "tumor type", "tumor region", "disease", "disease state",
                               "disease status"]:
                disease = value

            if key.lower() in ["background strain", "strain", "strain number", "mouse strain", "strain background",
                               "genetic background", "genotype", "genotype/variation",
                               "strain/background", "variation"]:
                genoType[key] = value

            if key.lower() in ["cell line", "cell", "cells pointed by barcodes", "cell line background",
                               "chicken line", "line"]:
                cellLine += value

            if key.lower() in ["cell_type", "cell-type", "cell type", "cell lineage"]:
                cellType += value

        organ_list = ['Bone', 'Cartilage', 'Fibrous joint', 'Cartilaginous joint', 'Synovial joint', 'Muscle', 'Tendon',
                      'Diaphragm', 'Heart', 'Bone marrow', 'Thymus', 'Spleen', 'Lymph node', 'CNS', 'Brain',
                      'Spinal cord', 'Ear', 'Eye', 'Integumentary system', 'Skin', 'Subcutaneous tissue', 'Breast',
                      'Mammary gland', 'Myeloid', 'Lymphoid', 'Nose', 'Nasopharynx', 'Larynx', 'Trachea', 'Bronchus',
                      'Lung', 'Mouth', 'Salivary', 'Tongue', 'Oropharynx', 'Laryngopharynx', 'Esophagus', 'Stomach',
                      'intestine', 'Appendix', 'Colon', 'Rectum', 'Anus', 'Liver', 'Biliary tract', 'Pancreas',
                      'Genitourinary', 'Kidney', 'Ureter', 'Bladder', 'Urethra', 'Uterus', 'Vagina', 'Vulva', 'Ovary',
                      'Placenta', 'Scrotum', 'Penis', 'Prostate', 'Testicle', 'Seminal vesicle', 'Pituitary', 'Pineal',
                      'Thyroid', 'Parathyroid', 'Adrenal', 'Islets of Langerhans']

        for organ in organ_list:
            if cellLine.lower().find(organ.lower()) != -1:
                sample.organ = organ
                break
            if cellType.lower().find(organ.lower()) != -1:
                sample.organ = organ
                break
            if tissue.lower().find(organ.lower()) != -1:
                sample.organ = organ
                break

        sample.antibody = antibody

        sample.treatment = treatment
        sample.disease = disease
        sample.cellLine = cellLine
        sample.cellType = cellType
        sample.genotype = genoType
        sample.tissue = tissue

        exclude = False
        for gse in sampleSeriesID:
            if gse in excludeGSEs:
                exclude = True

        if not exclude and sampleLibraryStrategy.lower() == 'chip-seq':
            results.append(sample)

            for gse_id in sample.series:
                GSE_GSM[gse_id].add(sample.id)

    queue.put((results, GSE_GSM))
    if db is not None:
        db.close()

settings = get_settings()
encode_pkl = settings['Encode']
roadmap_pkl = settings['Roadmap']
GGRmap_pkl = settings['GGR']

encodeGSE = load_obj(encode_pkl)

roadmapGSE = load_obj(roadmap_pkl)

GGRmapGSE = load_obj(GGRmap_pkl)

excludedGSE = set()

excludedGSE = excludedGSE.union(encodeGSE)
excludedGSE = excludedGSE.union(roadmapGSE)
excludedGSE = excludedGSE.union(GGRmapGSE)

db = sqlite3.connect('/home/tmhbxx3/scratch/GCF_OL/pkl/geoMetaData.db')
db.text_factory = str

query = db.execute('SELECT GSM_ID from GSM').fetchall()
localGSMs = list(set([x[0] for x in query]))

process = 20

queue = Queue()
chucksize = len(localGSMs)/(process-1) if process > 1 else len(localGSMs)
processes = []
for i in range(process):
    cur_geoGSMs = localGSMs[i*chucksize:(i+1)*chucksize]
    p = Process(target=chip_filter, args=(cur_geoGSMs, queue, excludedGSE))
    processes.append(p)
    p.start()

total_chipseqdatasets = []
GSE_GSM = {}

for i in range(process):
    results, cur_GSE_GSM = queue.get()
    total_chipseqdatasets += results
    GSE_GSM.update(cur_GSE_GSM)

for p in processes:
    p.join()

chip_db = sqlite3.connect('chipseq.db')
chip_db.execute('create table metadata(Data_ID text, Study_ID text, Title text, InstrumentModel text, RawData text, SequencingProtocol text, Species text, CellLine text, CellType text, Antibody text, Tissue text, Organ text, Characteristics text)')
chip_db.execute("CREATE INDEX index_metadata ON metadata (Data_ID);")

for sample in total_chipseqdatasets:
    chip_db.execute("insert into metadata values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                    (sample.id, ",".join(list(sample.series)), sample.title, sample.InstrumentID, sample.SRA,
                     sample.libraryStrategy, sample.organism, sample.cellLine, sample.cellType, json.dumps(sample.antibody),
                     sample.tissue, sample.organ, json.dumps(sample.characteristics)))
chip_db.commit()

chip_db.close()





#
# db = sqlite3.connect('/Users/boxia/Desktop/PycharmProjects/GEODataParser/GEOsearch/pkl/chipseq.db')
# db.text_factory = str
#
# query = db.execute('SELECT GSM_ID from GSM').fetchall()
# localGSMs = set([x[0] for x in query])
#
# process = 20
#
# queue = Queue()
# chucksize = len(localGSMs)/(process-1) if process > 1 else len(localGSMs)
# processes = []
#
# samples = {}
#
# for i in range(process):
#     cur_geoGSMs = localGSMs[i*chucksize:(i+1)*chucksize]
#     p = Process(target=feature_filter, args=(cur_geoGSMs, queue, features, features_begin, excludedGSM,
#                 type_seq, ignorecase, output_type,cwd,))
#     processes.append(p)
#     p.start()
#
# for i in range(process):
#     cur_samples = queue.get()
#     samples.update(cur_samples)
#
# for p in processes:
#     p.join()
#
# groupByGSE, excludedGSE, relatedSamples = SOFTQuickRelated(Human_Samples, output_type, type_seq,
#                                                                  GSEGSM_map, encode_remove, excludedGSE, cwd, process)