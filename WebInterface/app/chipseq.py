import sqlite3, json
from multiprocessing import Process, Queue
from sample import GSM
from input_search_utils import get_MetaInfo


def chip_filter(gsms, queue):
    db = sqlite3.connect('/Users/boxia/Desktop/PycharmProjects/GEODataParser/GEOsearch/pkl/geoMetaData.db')

    count = 0

    results = []

    for sampleName in gsms:
        sample = GSM(sampleName)

        sampleLibraryStrategy = ''

        info = get_MetaInfo(db, sample, count)

        count += 1

        for line in info:
            if line.startswith("!Sample_library_strategy"):
                sampleLibraryStrategy = line[line.find("=") + 1:].strip()
        if sampleLibraryStrategy.lower() == 'chip-seq':
            metadata = json.dumps(info)
            results.append((sampleName, metadata))

    queue.put(results)
    db.close()
    return


db = sqlite3.connect('/Users/boxia/Desktop/PycharmProjects/GEODataParser/GEOsearch/pkl/geoMetaData.db')
db.text_factory = str

query = db.execute('SELECT GSM_ID from GSM').fetchall()
localGSMs = list(set([x[0] for x in query]))

process = 20

queue = Queue()
chucksize = len(localGSMs)/(process-1) if process > 1 else len(localGSMs)
processes = []
for i in range(process):
    cur_geoGSMs = localGSMs[i*chucksize:(i+1)*chucksize]
    p = Process(target=chip_filter, args=(cur_geoGSMs, queue))
    processes.append(p)
    p.start()

total_chipseqdatasets = []

for i in range(process):
    total_chipseqdatasets += queue.get()

for p in processes:
    p.join()

chip_db = sqlite3.connect('chipseq.db')
chip_db.execute('create table GSM(GSM_ID text, metadata text)')
chip_db.execute("CREATE INDEX index_GSM ON GSM (GSM_ID);")

for gsm, metadata in total_chipseqdatasets:
    chip_db.execute("insert into GSM values(?, ?)", (gsm, metadata))
chip_db.commit()
chip_db.close()
