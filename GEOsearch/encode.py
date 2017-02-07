import pandas as pd

def encode_search(output_prefix, keywords, keywords_begin=(), type_seq='chip-seq', ignorecase=True, output_type='Homo sapiens'):
    seq_types = ['RNA-PET', 'ATAC-seq', 'genotyping by high throughput sequencing assay', 'HiC', 'RNA-seq',
                 'genetic modification followed by DNase-seq', 'ChIP-seq', '5C', 'siRNA knockdown followed by RNA-seq',
                 'eCLIP', 'TAB-seq', 'CRISPR genome editing followed by RNA-seq', 'RRBS', 'RIP-seq',
                 'shRNA knockdown followed by RNA-seq', 'iCLIP', 'RNA Bind-n-Seq',
                 'single cell isolation followed by RNA-seq', 'RAMPAGE', 'FAIRE-seq', 'DNase-seq',
                 'CAGE', 'PAS-seq', 'whole-genome shotgun bisulfite sequencing', "5' RLM RACE", 'Repli-seq',
                 'ChIA-PET', 'DNA-PET', 'microRNA-seq']
    seq_types_map = {}
    for type in seq_types:
        seq_types_map[type.lower()] = type

    if type_seq.lower() in seq_types_map:
        url = 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=fastq&assay_title='+\
              seq_types_map[type_seq.lower()]+'/metadata.tsv'
    else:
        url = 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=fastq/metadata.tsv'

    df = pd.read_csv(url, sep='\t')
    df = df[['File accession', 'Experiment accession', 'Controlled by', 'Assay', 'Biosample term name',
             'Biosample type', 'Biosample organism', 'File download URL', 'Platform', 'Experiment target', ]]

    df.columns = ['Sample_ID', 'Experiment_ID', 'Input', 'Type_Seq', 'Cell Line', 'Cell Type', 'Organism',
                  'Raw Data','Instrument_Model', 'Experiment target/antibody']

    case = False if ignorecase else True

    features = "|".join(keywords)
    df = df[df['Experiment target/antibody'].str.contains(features, case=case, na=False)]
    if ignorecase:
        features_begin = "|".join(keywords_begin).lower()
        df = df[df['Experiment target/antibody'].str.lower().str.startswith(features_begin, na=False)]
    else:
        features_begin = "|".join(keywords_begin)
        df = df[df['Experiment target/antibody'].str.startswith(features_begin, na=False)]

    df['Confidence'] = ['Very Confidence'] * len(df.index)
    df['Input_Description'] = ['indicated by encode'] * len(df.index)
    df[output_prefix.capitalize() + "_Description"] = df['Experiment target/antibody']

    samples_df = df.copy()
    samples_df = samples_df.set_index(['Sample_ID'])

    df = df[df['Organism'].str.contains(output_type, case=case, na=False)]
    df = df[df['Type_Seq'].str.contains(type_seq, case=case, na=False)]

    df = df.set_index(['Sample_ID'])

    return samples_df, df
