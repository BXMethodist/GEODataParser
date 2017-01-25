import pandas as pd
import pickle


def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


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
        try:
            id = id.strip()
            print id
            df = pd.read_csv("http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + id,
                             index_col=0)
            if table is None:
                table = df
            else:
                table = table.append(df)
        except:
            print id, " might not have SRA related information"
            failed.append(id)

    result_srr_gsm ={}
    result_srr_gse={}

    for srr in table.index.values:
        if srr in GSM_SRR_map:
            result_srr_gsm[srr] = GSM_SRR_map[srr]
            result_srr_gse[srr] = GSM_GSE_map[result_srr_gsm[srr]]
        else:
            print srr, " need to be updated"

    table['GSM_ID'] = pd.Series(result_srr_gsm)
    table['GSE_ID'] = pd.Series(result_srr_gse)

    table.to_csv(output_name, sep="\t")

    failed_file = open(output_name[:-4]+"_failed_id.txt", "w")

    for id in failed:
        failed_file.write(id+"\n")
    failed_file.close()










