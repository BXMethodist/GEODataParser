import pandas as pd
import pickle
from collections import defaultdict
import os

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


def GEO_query(names, output_name):
    # names could be a list or a file
    if isinstance(names, str) and os.path.exists(names) and os.path.isfile(names):
        list_names_obj = open(names, "r")
        list_names = [x.strip() for x in list_names_obj.readlines()]
        list_names_obj.close()
        names = list_names

    first_name = names[0]

    if first_name.startswith("GSE"):
        return query_GSE(names, output_name)
    else:
        return query_other(names, output_name)


def query_GSE(names, output_name):
    table = None
    result = defaultdict(set)
    GSM_GSE_map = load_obj("/home/tmhbxx3/scratch/XMLhttp/pickles/GSMGSE_map.pkl")

    GSMs = set()
    for name in names:
        GSMs = GSMs.union(GSM_GSE_map[name])
        for gsm in GSM_GSE_map[name]:
            result[gsm].add(name)

    for gsm in GSMs:
        df = pd.read_csv("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?sp=runinfo&acc=" + gsm + "&retmode=txt",
                         index_col=29)
        if table is None:
            table = df
        else:
            table = table.append(df)
    table['GSE_ID'] = pd.Series(result)
    table.to_csv(output_name + ".txt", sep="\t")
    return table

def query_other(names, output_name):
    table = None
    for name in names:
        df = pd.read_csv("https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?sp=runinfo&acc="+name+"&retmode=txt", index_col=29)

        if table is None:
            table = df
        else:
            table = table.append(df)

    GSM_GSE_map = load_obj("/home/tmhbxx3/scratch/XMLhttp/pickles/GSMGSE_map.pkl")

    result = {}
    GSMs = table.index.values
    for gsm in GSMs:
        result[gsm] = GSM_GSE_map[gsm]
    table['GSE_ID'] = pd.Series(result)
    table.to_csv(output_name+".txt", sep="\t")
    return table





