#!/usr/bin/env python3
import os, re
import pandas as pd

#script adding column with whole description from batchentrez database to dataframe with tAI, protein ID and protein fammily name

batchentrez = "/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/batchentrez/"
protein_ids_folder = '/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/protein_ids'
diff_exp_proteins = '/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/diff_exp_proteins_tAI.txt'
output_path = '/home/norbert_s/tAI/all_chunks/bialka_dla_Magdy/final_results.tsv'
df = pd.read_csv(diff_exp_proteins, sep='\t')


leng = 0
file_count = 0
protein_description_df = pd.DataFrame(columns=['protein_id', 'whole_description'])
for filename in os.listdir(batchentrez):
    protein_patterns = []
    filename_shorted = filename.split('.')[0] + '.modified'
    file_count += 1
    with open(os.path.join(batchentrez, filename), "r") as f:
        f = f.read()
        f = f.strip().split("\n")

        description = " ".join(f)
        all_descriptions = re.split(r"\s\d+\.\s", description)
        all_descriptions = [description.strip() for description in all_descriptions]
        all_descriptions = [re.sub(r'^\d+\.\s*', '', item) for item in all_descriptions]
        # print(all_descriptions)
        for idx, id in enumerate(df.protein_id):
            matching_desc = [desc for desc in all_descriptions if id in desc]
            if matching_desc:
                df.loc[idx, 'whole_description'] = matching_desc[0]
        # [print(all_descriptions)]
        leng =+ len(all_descriptions)
        # print(all_descriptions)
print('przed usunieciem NAN\n')
print(df.apply(len))
# print(df.iloc[0:10, 2])
# print(df.iloc[0:10, 3])
df = df.dropna()
print('po usunieciu NAN')
print(df.apply(len))


df.to_csv(output_path, sep='\t', index=False)


