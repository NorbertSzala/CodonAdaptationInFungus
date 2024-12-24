#!/bin/bash

# function to check for errors after each command
check_error() {
    if [ $? -ne 0 ]; then
        echo "Error in $1"
        exit 1
    fi
}


#setting temporary folder
export TMPDIR=/home/norbert_s/tmp

#setting pathways to genomic and cds seqs
genomic_folder_path="/db/183_fungi_mai_2022/genomic_fna"
cds_folder_path="/db/183_fungi_mai_2022/cds"

# genomic_folder_path="/home/norbert_s/genome"
# cds_folder_path="/home/norbert_s/183_CDS"

#checking if input folders exist
if [ ! -d "$genomic_folder_path"  ]; then
    echo "$genomic_folder_path doesn't exist!"
    exit 1
fi

if [ ! -d "$cds_folder_path"  ]; then
    echo "$cds_folder_path doesn't exist!"
    exit 1
fi



# Creating arrays with genomic and CDS sequence files
genomic_seq_files=($(find "$genomic_folder_path" -type f -name "*.fna"))
cds_seq_files=($(find "$cds_folder_path" -type f -name "*.fna"))

tRNAscan-SE_files=($(find "$tRNAscan-SE_output" -type f -name "*.fna")) # to ma isc do petli ktora iteruje po wszystkich folderach. Dodaj prefix do koncowki folderu


#creating table to store files with prefix
declare -A genomic_files_map

#initialising map 
echo "Creating map started"
for genomic_seq in "${genomic_seq_files[@]}"; do
    genomic_prefix=$(basename "$genomic_seq" | cut -c1-13) #creating prefix whis is common with assemblyID
    genomic_files_map["$genomic_prefix"]+="$genomic_seq"
done
echo "Creating map done"    


file_pairs=()
echo "Creating pairs started"
for cds_seq in "${cds_seq_files[@]}"; do
    cds_prefix=$(basename "$cds_seq" | cut -c1-13) #creating prefix whis is common with assemblyID
        if [[ -n "${genomic_files_map[$cds_prefix]}" ]]; then
        for genomic_seq in ${genomic_files_map[$cds_prefix]}; do #combining genomic seq and cds seq to pairs according to common prefix (assemblyID)
            file_pairs+=("$genomic_seq $cds_seq")
        done
    fi
done
echo "Creating pairs done"
echo "pairs: $file_pairs"   




# process_pair() {
#     pair=$1
#     counter=$2
#     output_folder=$3

#     echo "Processing pair number $counter: $pair"

#     # extracting names of files from pair
#     genomic_seq=$(echo "$pair" | cut -d ' ' -f 1)
#     cds_seq=$(echo "$pair" | cut -d ' ' -f 2)

#     # creating full path to files
#     genomic_seq_full_path="$genomic_folder_path/$(basename "$genomic_seq")"
#     cds_seq_full_path="$cds_folder_path/$(basename "$cds_seq")"

#     prefix_inside="$(basename "$genomic_seq" | cut -c1-13)"
#     echo "hej"

#     # Main operations
#     # # # # Process CDS sequence using other tools sequentially
#     echo "Python script count_genes_coding_tRNA.py started for $pair"
#     python3 count_genes_coding_tRNA.py tRNAscan-SE_output_$prefix_inside counted_genes_coding_tRNA_$prefix_inside.txt
#     check_error "count_genes_coding_tRNA.py for $cds_seq"


#     # echo "CodonM started for $pair"
#     # perl codonM "$cds_seq_full_path" codonM_output_$prefix_inside
#     # check_error "codonM for $cds_seq"

#     # echo "CodonW started for $pair"
#     # codonw "$cds_seq_full_path" codonw_output_$prefix_inside codons_NC_bulk_output_$prefix_inside -nomenu -nowarn -silent -enc -gc -gc3s -L_aa
#     # check_error "codonw for $cds_seq"


#     # echo "Rscript started for $pair"
#     # echo "output folder: $output_folder"
#     # Rscript tAI.R counted_genes_coding_tRNA_$prefix_inside.txt codonM_output_$prefix_inside codonw_output_$prefix_inside $output_folder/each_protein_tAI_$prefix_inside.txt $output_folder/summary_tAI_$prefix_inside.txt
#     # check_error "Rscript tAI.R for $cds_seq"

#     # # # Clean up intermediate files
#     # rm -f counted_genes_coding_tRNA_$prefix_inside.txt codonM_output_$prefix_inside codonw_output_$prefix_inside codons_NC_bulk_output_$prefix_inside

#     echo "Processing for $genomic_seq and $cds_seq completed. It was $counter th pair"
#     echo " "
#     echo " "
# }




base_folder="/home/norbert_s/tAI/lazy_output/output_shorted_pipeline_tAI"
for folder_index in {1..10}; do
    tRNAscan-SE_folder="${base_folder}${folder_index}/tRNAscan-SE"
    counted_genes_coding_tRNA_folder="${base_folder}${folder_index}/counted_genes_coding_tRNA_folder"
    mkdir -p "$counted_genes_coding_tRNA_folder"

    if [ -d "$tRNAscan-SE_folder" ]; then
        for file_path in "$tRNAscan-SE_folder"/*; do
            if [ -f "$file_path" ]; then
                prefix=$(basename "$file_path" | tail -c 13)

                echo "Python script count_genes_coding_tRNA.py started for $file_path"
                python3 /home/norbert_s/tAI/count_genes_coding_tRNA.py "$file_path/tRNAscan-SE_output_$prefix" "$counted_genes_coding_tRNA_folder/counted_genes_coding_tRNA_$prefix.txt"
                check_error "count_genes_coding_tRNA.py for $prefix"
                echo " "
                echo " "
            fi
        done
    fi
done

