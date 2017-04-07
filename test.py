import os
import pandas as pd

old_df = pd.read_excel('./GEOsearch/search_output/paper/Table_S1_ChipSeqPair_H3K4me3_SearchResults.xlsx', index_col=0)

new_df = pd.read_csv('./GEOsearch/search_output/Search_ResultsampleWithH3K4me3.csv', index_col=0)

old_df['Organ'] = new_df.ix[old_df.index, 'Organ']
old_df['Tissue'] = new_df.ix[old_df.index, 'Tissue']

old_df.to_excel('Table_S1_ChipSeqPair_H3K4me3_SearchResults.xlsx')





