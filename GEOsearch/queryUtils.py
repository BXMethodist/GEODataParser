import pandas as pd
import pickle


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


def encode_metadata(id):
    try:
        url = 'https://www.encodeproject.org/metadata/type=Experiment&files.accession='+id+'/metadata.tsv'
        df = pd.read_csv(url, sep='\t')
        df = df[['File accession', 'Experiment accession', 'Read length', 'Run type', 'Paired with', 'File download URL']]
        df.columns = ['Run_ID', 'Experiment_ID', 'Read length', 'Run type', 'Paired with', 'File download URL']
        df = df[df['Run_ID'] == id]
        df = df.set_index(['Run_ID'])
        return df
    except:
        return None


def geo_metadata(id):
    try:
        url = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + id
        df = pd.read_csv(url)
        df = df[['Run', 'avgLength', 'LibraryLayout', 'download_path']]
        df.columns = ['Run_ID', "Read length", "Run type", 'File download URL']
        df = df.set_index(['Run_ID'])
        return df
    except:
        return None


def GEO_query(names, output_name, GSM_GSE_pkl, GSM_SRR_pkl):
    # names could be a list

    GSM_GSE_map = load_obj(GSM_GSE_pkl)
    GSM_SRR_map = load_obj(GSM_SRR_pkl)

    query_ids = set()

    for name in names:
        if name.startswith("GSE"):
            for gsm in GSM_GSE_map[name]:
                query_ids.add(gsm)
        else:
            query_ids.add(name)

    table = None

    failed = []

    for id in query_ids:
        id = id.strip()
        if id.startswith('ENC'):
            df = encode_metadata(id)
            if df is not None and len(df.index) > 0:
                if table is None:
                    table = df
                else:
                    table = table.append(df)
        else:
            df = geo_metadata(id)
            if df is not None and len(df.index) > 0:
                if table is None:
                    table = df
                else:
                    table = table.append(df)
        if df is None or len(df.index) == 0:
            print id, " might not have SRA related information"
            failed.append(id)

    result_srr_gsm = {}
    result_srr_gse = {}

    if table is not None:
        for srr in table.index.values:
            if srr in GSM_SRR_map:
                result_srr_gsm[srr] = GSM_SRR_map[srr]
                result_srr_gse[srr] = ','.join(list(GSM_GSE_map[result_srr_gsm[srr]]))
            else:
                print srr, " need to be updated"

        table['GSM_ID'] = pd.Series(result_srr_gsm)
        table['GSE_ID'] = pd.Series(result_srr_gse)

        table.to_csv(output_name, sep="\t")

    failed_file = open(output_name[:-4]+"_failed_id.txt", "w")

    for id in failed:
        failed_file.write(id+"\n")
    failed_file.close()









