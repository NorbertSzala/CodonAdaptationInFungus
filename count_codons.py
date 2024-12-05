#!/usr/bin/python3

'''
Program reads table from tRNAscan-SE output and count tRNA genes that code codons. Codons in TCAG order. Output is in format needed in gtAI (R) or cubar (R package).
I used reversed translation, because in tRNAscan-SE output is given sequence of Anticodon
Output is number of genes that code specific tRNA which have determined anticodon.
Note: input Anti Codon is in DNA (as seq in input tRNAscan-SE)
'''

from collections import OrderedDict

# input_path = '/home/norbert/IBB/skrypty/tRNAscan-SE_postia_output'
input_path = '/home/norbert/IBB/skrypty/tRNAscan-SE_diff_output'
# #input_path = '/home/norbert_s/IBB/skrypty/tRNAscan-SE_postia_output'
#
output_path = '/home/norbert/IBB/skrypty/count_codons.txt'
# #output_path = '/home/norbert_s/IBB/skrypty/count_codons.txt'

# input_path = '/home/norbert/IBB/skrypty/example_seq.txt'

translate = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C'} #dict to translate nts.
codons_order = '''
TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG
CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG
ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG
GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG
'''.split()

def read_anticodons(path: str) -> list:
    """Reads anticodons from tRNAscan-SE output file."""
    try:
        with open(path, 'r') as file:
            lines = file.readlines()
        # Extract anticodons from the 6th column, skipping headers

            return [line.strip().upper().split('\t')[5] for line in lines[3:]]

    except Exception as e:
        print(f"Error reading file: {e}")
        return []

def rev_tran_anticodon(anticodon: str) -> str:
    """Translates an anticodon to a codon (reverse complement)."""
    try:
        return ''.join(translate[nt] for nt in anticodon)[::-1]
    except KeyError:
        return None  # Invalid anticodon


def count_codons(anticodons: list, only_count_anticodon: bool = True) -> OrderedDict:
    """Counts codons in the given anticodon list."""
    codons_count = {}

    for anticodon in anticodons:
        print(anticodon)
        if only_count_anticodon:
            codons_count[anticodon] = codons_count.get(anticodon,0) + 1  # check if anticodon is present in dict. If not, it adds 1 value

        else:
            codon = rev_tran_anticodon(anticodon)
            if codon:
                codons_count[codon] = codons_count.get(codon, 0) + 1 #check if codon is present in dict. If not, it adds 1 value
        print(sum(codons_count.values()))

    # Ensure all 64 codons are present. If not, lacked codons value is set as 0.
    for codon in codons_order:
        codons_count.setdefault(codon, 0)

    # Return the OrderedDict and invalid nucleotides as a tuple
    return OrderedDict((codon, codons_count[codon]) for codon in codons_order) #return sorted dick in TCAG order

def write_output(codons_count, path):
    """Writes the codon counts to the output file."""
    try:
        with open(path, 'w') as file:
            file.writelines([f"{count}\n" for count in codons_count.values()])
    except Exception as e:
        print(f"Error writing file: {e}")

def main():
    anticodons = read_anticodons(input_path)
    if anticodons:  # Ensure anticodons were successfully read
        codons_count = count_codons(anticodons, False)

        write_output(codons_count, output_path)
        print(f"Codon counts written to {output_path}")
    else:
        print("No valid anticodons found.")

if __name__ == "__main__":
    main()
