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
# genomic_folder_path="/db/183_fungi_mai_2022/genomic_fna"
# cds_folder_path="/db/183_fungi_mai_2022/cds"

genomic_folder_path="/home/norbert_s/genome"
cds_folder_path="/home/norbert_s/183_CDS"

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


#making ten folders where output data will be saved. It is easier to quality controll and debugging
for i in {1..10}; do
    mkdir "output_shorted_pipeline_tAI$i"
done


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


process_pair() {
    pair=$1
    counter=$2
    output_folder=$3

    echo "Processing pair number $counter: $pair"

    # extracting names of files from pair
    genomic_seq=$(echo "$pair" | cut -d ' ' -f 1)
    cds_seq=$(echo "$pair" | cut -d ' ' -f 2)

    # creating full path to files
    genomic_seq_full_path="$genomic_folder_path/$(basename "$genomic_seq")"
    cds_seq_full_path="$cds_folder_path/$(basename "$cds_seq")"

    prefix_inside="$(basename "$genomic_seq" | cut -c1-13)"

    # Main operations
    start_time=$(date +%s)
    echo "tRNAscan-SE started for $pair"
    tRNAscan-SE -o tRNAscan-SE_output_$prefix_inside "$genomic_seq_full_path"
    # tRNAscan_pid=$!  # Store the process ID of tRNAscan-SE

    # Wait for tRNAscan-SE to finish
    # wait $tRNAscan_pid
    check_error "tRNAscan-SE for $genomic_seq"

    echo "tRNAscan-SE done"
    end_time=$(date +%s)
    echo "tRNAscan-SE execution time: $((end_time - start_time)) seconds"



    # # # Process CDS sequence using other tools sequentially
    echo "Python script count_genes_coding_tRNA.py started for $pair"
    python3 count_genes_coding_tRNA.py tRNAscan-SE_output_$prefix_inside counted_genes_coding_tRNA_$prefix_inside.txt
    check_error "count_genes_coding_tRNA.py for $cds_seq"

    # rm tRNAscan-SE_output_$prefix_inside

    echo "CodonM started for $pair"
    perl codonM "$cds_seq_full_path" codonM_output_$prefix_inside
    check_error "codonM for $cds_seq"

    echo "CodonW started for $pair"
    codonw "$cds_seq_full_path" codonw_output_$prefix_inside codons_NC_bulk_output_$prefix_inside -nomenu -nowarn -silent -enc -gc -gc3s -L_aa
    check_error "codonw for $cds_seq"


    echo "Rscript started for $pair"
    echo "output folder: $output_folder"
    Rscript tAI.R counted_genes_coding_tRNA_$prefix_inside.txt codonM_output_$prefix_inside codonw_output_$prefix_inside $output_folder/each_protein_tAI_$prefix_inside.txt $output_folder/summary_tAI_$prefix_inside.txt
    check_error "Rscript tAI.R for $cds_seq"

    # Clean up intermediate files
    rm -f counted_genes_coding_tRNA_$prefix_inside.txt codonM_output_$prefix_inside codonw_output_$prefix_inside codons_NC_bulk_output_$prefix_inside

    echo "Processing for $genomic_seq and $cds_seq completed. It was $counter th pair"
    echo " "
    echo " "
}

#Eexport variables to parallel function
export -f process_pair
export -f check_error
export genomic_folder_path
export cds_folder_path
export files_per_folder=19
export output_base_folder="output_shorted_pipeline_tAI"
num_pairs=${#file_pairs[@]}


# Execute the process in parallel on 80 cores
echo "Processing started"

counter=1
for pair in "${file_pairs[@]}"; do
    folder_index=$((counter / files_per_folder + 1))
    output_folder="$output_base_folder$folder_index"
    # Przekazywanie argumentów do parallel
    echo "$pair $counter $output_folder"
    counter=$((counter + 1))
done | parallel -j 80 'process_pair {1} {2} {3}'



# DO TEJ PORY DZIALA LECZ BEZ DODANEJ OBSLUGI 80 RDZENIOW
# #creating table with file pairs
# file_pairs=("${file_pairs[@]}")

# output_folders=("output_shorted_pipeline_tAI1" "output_shorted_pipeline_tAI2" "output_shorted_pipeline_tAI3" "output_shorted_pipeline_tAI4" "output_shorted_pipeline_tAI5" 
#                 "output_shorted_pipeline_tAI6" "output_shorted_pipeline_tAI7" "output_shorted_pipeline_tAI8" "output_shorted_pipeline_tAI9" "output_shorted_pipeline_tAI10")

# files_per_folder=19
# num_folders=${#output_folders[@]}  # Number of folders in the output folder
# num_pairs=${#file_pairs[@]}  # Number of file pairs

# counter=1

# for pair in "${file_pairs[@]}"; do

#     # Calculate folder index based on how many pairs have been processed
#     folder_index=$((counter / files_per_folder))
#     # Check if we have more pairs than available folders
#     if ((folder_index >= num_folders)); then
#         break
#     fi

#     # Path to the output folder
#     output_folder=${output_folders[$folder_index]}  # Folder index corresponds to the right folder
#     output_file="$output_folder/output_$counter.txt"

#     #extracting names of files from pair
#     genomic_seq=$(echo "$pair" | cut -d ' ' -f 1)
#     cds_seq=$(echo "$pair" | cut -d ' ' -f 2)

#     #creating full path to files
#     genomic_seq_full_path="$genomic_folder_path/$(basename "$genomic_seq")"
#     cds_seq_full_path="$cds_folder_path/$(basename "$cds_seq")"



#     # Main operations
#     tRNAscan-SE -o tRNAscan-SE_output_$genomic_prefix "$genomic_seq_full_path" &
#     tRNAscan_pid=$!  # Store the process ID of tRNAscan-SE

#     # Wait for tRNAscan-SE to finish
#     wait $tRNAscan_pid
#     check_error "tRNAscan-SE for $genomic_seq"

#     # Process CDS sequence using other tools sequentially
#     python3 count_genes_coding_tRNA.py tRNAscan-SE_output_$genomic_prefix counted_genes_coding_tRNA_$genomic_prefix.txt
#     check_error "count_genes_coding_tRNA.py for $cds_seq"

#     rm tRNAscan-SE_output_$genomic_prefix

#     perl codonM "$cds_seq_full_path" codonM_output_$genomic_prefix
#     check_error "codonM for $cds_seq"

#     codonw "$cds_seq_full_path" codonw_output_$genomic_prefix codons_NC_bulk_output_$genomic_prefix -nomenu -nowarn -silent -enc -gc -gc3s -L_aa
#     check_error "codonw for $cds_seq"

#     Rscript tAI.R counted_genes_coding_tRNA_$genomic_prefix.txt codonM_output_$genomic_prefix codonw_output_$genomic_prefix each_protein_tAI_$genomic_prefix.txt summary_tAI_$genomic_prefix.txt
#     check_error "Rscript tAI.R for $cds_seq"

#     # Clean up intermediate files
#     rm -f counted_genes_coding_tRNA_$genomic_prefix.txt codonM_output_$genomic_prefix codonw_output_$genomic_prefix codons_NC_bulk_output_$genomic_prefix

#     echo "Processing for $genomic_seq and $cds_seq completed"
    
#     # Zwiększamy licznik
#     counter=$((counter + 1))
# done

# echo "Processing completed."
