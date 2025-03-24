#!/usr/bin/env python3

import pandas as pd
'''
concatenate df with average tAI for each 183 genomes and df with assembly id 
'''

tAI_df_path = '/home/norbert_s/tAI/all_chunks/concatenated_summary_protein_file_shorted.txt'
main_df_path = '/home/norbert_s/tAI/all_chunks/assembly_id_183.txt'

tAI_df = pd.read_csv(tAI_df_path, delim_whitespace=True)
main_df = pd.read_csv(main_df_path, header = None)

main_df['assembly_ID'] = main_df[0].str.split('.').str[0]
main_df = main_df[['assembly_ID']]


merged_df = pd.merge(main_df, tAI_df, on='assembly_ID',  how = 'left')

merged_df = merged_df.set_index('assembly_ID').reindex(main_df['assembly_ID']).reset_index()
merged_df.to_csv('/home/norbert_s/tAI/all_chunks/assembly_id_tAI_183.txt', sep='\t', index=False)

