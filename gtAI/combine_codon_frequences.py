#!/usr/bin/env python3
#Short script to count total number of different codons in table
import pandas as pd
import os


# path = '/home/norbert/IBB/skrypty/codon_frequencies_postia'
path = '/home/norbert_s/gtAI/codon_frequencies_postia'

# path_to_save = '/home/norbert/IBB/skrypty/combined_codons_frequencies.txt'
path_to_save = '/home/norbert_s/gtAI/combined_codons_frequencies.txt'

def combine_codons_frequencies(path, path_to_save):

    if not os.path.exists(path):
        print(f'Error: The file {path} does not exist!')
        return

    data_temp = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line = list(map(int, line))
            data_temp.append(line) #firstly create one large list of lists

        df = pd.DataFrame(data_temp)

        total_codons = df.sum(axis=0) #sum values from each column to count total amount of codons
        total_codons = pd.DataFrame(total_codons) #save values as dataframe
        total_codons = total_codons.T #transpose a df
        total_codons.to_csv(path_to_save, index=False, header=False, sep='\t') #save as csv (tsv) file

combine_codons_frequencies(path, path_to_save)
