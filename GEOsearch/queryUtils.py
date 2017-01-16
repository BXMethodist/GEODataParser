import pandas as pd
import pickle
from collections import defaultdict
import os

def load_obj(name):
    with open(name, 'rb') as f:
        return pickle.load(f)


def GEO_query(names, output_name):
    # names could be a list
    GSEs = []

    non_GSEs = []
    for name in names:
        if name.startswith("GSE"):
            GSEs.append(name)
        else:
            non_GSEs.append(name)

    table1 = query_GSE(GSEs)
    table2 = query_other(non_GSEs)

    if table1 is None:
        table2.to_csv(output_name, sep="\t")
        return table2
    elif table2 is None:
        table1.to_csv(output_name, sep="\t")
        return table1
    else:
        table = table1.append(table2)
        table.to_csv(output_name, sep="\t")
        return table


def query_GSE(names):
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
    return table

def query_other(names):
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
    return table





