# tRNAscan-SE programme results:

## command:
tRNAscan-SE -H -o tRNAscan-SE_postia_output -f tRNAscan-SE_postia_structure -s tRNAscan-SE_postia_isospecific -m tRNAscan-SE_postia_stats -j tRNAscan-SE_postia_gff -a tRNAscan-SE_postia_seqs /db/183_fungi_mai_2022/genomic_fna/GCA_000006255.1_Postia_placenta_V1.0_genomic.fna

## Output files:

### tRNAscan-SE_postia_output
main file with holistic results - table with seqs names and predicted tRNA with anti codon and nucleotide index where tRNA begins and ends. It consist aldo some score indicators with adnotation of being pseudogene


### tRNAscan-SE_postia_structure
Informations about structure in sequence. It consist of anticodon, score, type (aa) or 'graphical' scheme.


### tRNAscan-SE_postia_isospecific
Informs about what other aminoacids (except for special one preferred aa) can be transported by tRNA. Other words - aminoacid preferencess in tRNA


### tRNAscan-SE_postia_gff
Classic gff format table with informations from tRNAscan-SE_postia_output


### tRNAscan-SE_postia_seqs
Fasta format. Sequences of tRNA with index and predicted transported aminoacid, codon/anticodon i dont know, score and length. Really simply form. 


### tRNAscan-SE_postia_stats
Sequences read:         1242
Seqs w/at least 1 hit:  147
Bases read:             90891856 (x2 for both strands)
Bases in tRNAs:         25480
tRNAs predicted:        296
Av. tRNA length:        86
Predicted pseudogenes:	28

