#!/usr/bin/env python3
#this script makes two plots (counts and percentage) of aminoacids to compare data with general amino acid composihttps://www.ebi.ac.uk/uniprot/TrEMBLstats
#as input 1 (path) is needed a ?table? with frequency, names and counts of codons calculated by codonw.

import matplotlib.pyplot as plt
import numpy as np
# path = '/home/norbert/IBB/skrypty/codons_NC_bulk_totals'
# path_to_plots = '/home/norbert/IBB/skrypty/'
path = '/home/norbert_s/gtAI/codons_NC_bulk_totals'
path_to_plots = '/home/norbert_s/gtAI/'

def aminoacidslayout(path):
    codons_count = {}
    aminoacids_count = {}

    codon_table = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}
    codon_table = dict(sorted(codon_table.items(), key = lambda item: item[0])) #sort in alphabetical order. Two dicts have to be in the same order

    with open(path, 'r') as f:
        words = [word for word in f.read().strip().split()]

        for word in words:
            if word[:3].isalpha() and word[3:].isdigit() and len(word) > 3: #Extracting codon sequences and occurrence counts
                codon = word[:3]
                count = int(word[3:])
                if codon not in codons_count:
                    codons_count[codon] = count
                else:
                    codons_count[codon] += count

        codons_count = dict(sorted(codons_count.items(), key = lambda item: item[0])) #another dict sorting

        if len(codon_table) == len(codons_count): #checking if dicts have the same amount of values
            for codon_intable, codon_incount in zip(codon_table.keys(), codons_count.keys()):
                if codon_intable == codon_incount:
                    amino_acid = codon_table[codon_intable]
                    if amino_acid not in aminoacids_count:
                        aminoacids_count[amino_acid] = codons_count[codon_incount] #creating new dicti with name of AA and occurrence count
                    else:
                        aminoacids_count[amino_acid] += codons_count[codon_incount]

        aminoacids_count = dict(sorted(aminoacids_count.items(), key = lambda item: item[1], reverse = True)) #sorting ascending

        total = sum([value for value in aminoacids_count.values()])
        aminoacids = [aa for aa in aminoacids_count.keys()]
        counts = [count for count in aminoacids_count.values()]
        percentage = [round(count/total*100, 2) for count in aminoacids_count.values()]


        plt.bar(aminoacids, counts, width = 0.5)
        plt.xticks(aminoacids, rotation = 90)
        plt.xlabel('Amino Acids')
        plt.ylabel('Counts')
        plt.title('Frequency of occuring amino acids')
        plt.tight_layout()
        plt.savefig(path_to_plots+"amino_acids_counts.png", format="png", dpi=300, bbox_inches="tight")
        plt.show()


        plt.bar(aminoacids, percentage, width = 0.5)
        plt.xticks(aminoacids, rotation = 90)
        plt.xlabel('Amino Acids')
        plt.ylabel('Percentage')
        plt.title('Percentage frequency of occuring amino acids')
        plt.savefig(path_to_plots+"amino_acids_percentage.png", format="png", dpi=300, bbox_inches="tight")
        plt.show()



aminoacidslayout(path)

