#!/usr/bin/python3

import pandas as pd

# Paths
protein_tAI = "/home/norbert_s/tAI/all_chunks/concatenated_each_protein_file_shorted.txt"
pointed_proteins = "/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/all_proteins.txt"
output = "/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/diff_exp_proteins_tAI.txt"

# Reading data
df_tAI = pd.read_csv(protein_tAI, sep="\t")
df_proteins = pd.read_csv(pointed_proteins, sep=" ", header=None)
# print(df_tAI.iloc[:, 1])
df_proteins.columns = ['protein_id', 'feature']

# Sorting data frames by protein_id column
df_tAI_sorted = df_tAI.sort_values(by="protein_id")
df_proteins_sorted = df_proteins.sort_values(by="protein_id")

df_proteins_sorted = df_proteins_sorted.dropna(subset=['protein_id'])
df_tAI = df_tAI.dropna(subset=['protein_id'])

print(df_proteins.head())
print(df_proteins.shape)
print(df_tAI.head())
print(df_tAI.shape)

# Merging data by protein_id column
merged_df = pd.merge(df_tAI_sorted, df_proteins_sorted, on="protein_id", how="inner")

merged_selected = merged_df.iloc[:, [0,1,2,3]]
print(merged_selected.head())
merged_selected.to_csv(output, index=False, sep="\t")
