import pandas as pd

# Ścieżki do plików
protein_tAI = "/home/norbert_s/tAI/all_chunks/concatenated_each_protein_file_shorted.txt"
pointed_proteins = "/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/all_proteins.txt"
output = "/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/diff_exp_proteins_tAI.txt"

# Wczytanie danych do pandas DataFrame
df_tAI = pd.read_csv(protein_tAI, sep="\t")  # Zakładając, że plik jest oddzielony tabulatorami
print(df_tAI.head())
df_proteins = pd.read_csv(pointed_proteins, sep="\t", header=None)

df_proteins.columns = ['protein_id']

#
# # Posortowanie obu DataFrame'ów po kolumnie protein_id
df_tAI_sorted = df_tAI.sort_values(by="protein_id")
df_proteins_sorted = df_proteins.sort_values(by="protein_id")
print(df_tAI.head())
print(df_proteins.head())
# # Połączenie danych na podstawie kolumny protein_id
merged_df = pd.merge(df_tAI_sorted, df_proteins_sorted, on="protein_id", how="inner")
merged_selected = merged_df.iloc[:, [1, 2]]

# # Zapisanie wynikowego DataFrame do pliku
merged_selected.to_csv(output, index=False, sep="\t")
