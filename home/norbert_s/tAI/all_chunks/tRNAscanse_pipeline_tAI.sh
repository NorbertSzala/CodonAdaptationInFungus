#!/bin/bash

### Script Summary:
# This script automates the process of running tRNAscan-SE on pairs of genomic and CDS (coding sequence) files,
# organizing the results into multiple output folders. Main goal is to count tRNA in cds. The script is designed to handle large datasets by
# dividing the processing into smaller chunks, with each chunk processed in parallel to optimize execution time.
# This is the first of two steps to calculate tAI. The process has been divided into two smaller parts
# to optimize error handling and execution time, as tRNAscan-SE can be time-consuming to run.

#Second and last part of script is /home/norbert_s/tAI/all_chunks/other_analysis_pipeline_tAI.sh

# 1. Sets up paths and verifies the existence of the genomic and CDS sequence folders.
# 2. Creates 10 output folders to store the results of tRNAscan-SE.
# 3. Maps genomic files by prefixes (assembly IDs) and creates pairs with CDS sequences.
# 4. Processes each genomic and CDS pair by running tRNAscan-SE and storing results in corresponding folders.
# 5. Uses parallel processing for efficient handling of file pairs, processing up to 19 pairs per folder.
# 6. Includes error handling to stop execution if any critical step fails.

### Key Variables:
# genomic_folder_path - Path to the folder containing genomic sequence files.
# cds_folder_path - Path to the folder containing CDS sequence files.
# output_folder - Folder where results from tRNAscan-SE will be stored.
# file_pairs - Array holding pairs of genomic and CDS files.
# files_per_folder - The number of file pairs processed per folder (set to 19 in this script).
# tRNAscan-SE - Tool used for identifying tRNA genes in genomic sequences.

### Example Usage:
# To use this script, ensure that the paths to the genomic and CDS sequence folders are correct,
# and that tRNAscan-SE is installed and accessible. The script automatically divides the workload into chunks,
# processes them, and organizes the results into corresponding folders.
# Function to check for errors after each command


check_error() {
    if [ $? -ne 0 ]; then
        echo "Error in $1"
        exit 1
    fi
}

# Setting temporary folder
export TMPDIR=/home/norbert_s/tmp

# Setting pathways to genomic and cds seqs
genomic_folder_path="/db/183_fungi_mai_2022/genomic_fna"
cds_folder_path="/home/norbert_s/tAI/all_chunks/cds"

# Checking if input folders exist
if [ ! -d "$genomic_folder_path" ]; then
    echo "$genomic_folder_path doesn't exist!"
    exit 1
fi

if [ ! -d "$cds_folder_path" ]; then
    echo "$cds_folder_path doesn't exist!"
    exit 1
fi

# Creating arrays with genomic and CDS sequence files
genomic_seq_files=($(find "$genomic_folder_path" -type f -name "*.fna"))
cds_seq_files=($(find "$cds_folder_path" -type f -name "*.fna"))

# Making ten folders where output data will be saved
for i in {1..10}; do
    mkdir -p "output_shorted_pipeline_tAI$i"
done

# Creating a table to store files with the prefix
declare -A genomic_files_map

# Initializing map
echo "Creating map started"
for genomic_seq in "${genomic_seq_files[@]}"; do
    genomic_prefix=$(basename "$genomic_seq" | cut -c1-13) # Creating prefix which is common with assemblyID
    genomic_files_map["$genomic_prefix"]+="$genomic_seq"
done
echo "Creating map done"

file_pairs=()
echo "Creating pairs started"
for cds_seq in "${cds_seq_files[@]}"; do
    cds_prefix=$(basename "$cds_seq" | cut -c1-13) # Creating prefix which is common with assemblyID
    if [[ -n "${genomic_files_map[$cds_prefix]}" ]]; then
        for genomic_seq in ${genomic_files_map[$cds_prefix]}; do # Combining genomic seq and cds seq to pairs according to common prefix (assemblyID)
            file_pairs+=("$genomic_seq $cds_seq")
        done
    fi
done
echo "Creating pairs done"
echo "pairs: $file_pairs"

process_pair() {
    pair=$1
    counter=$2
    output_folder=$3

    echo "Processing pair number $counter: $pair"

    # Extracting names of files from pair
    genomic_seq=$(echo "$pair" | cut -d ' ' -f 1)
    cds_seq=$(echo "$pair" | cut -d ' ' -f 2)

    # Creating full path to files
    genomic_seq_full_path="$genomic_folder_path/$(basename "$genomic_seq")"
    cds_seq_full_path="$cds_folder_path/$(basename "$cds_seq")"

    prefix_inside="$(basename "$genomic_seq" | cut -c1-13)"

    # Main operations (for demo, echo instead of tRNAscan-SE)
    echo "tRNAscan-SE started for $pair"
    tRNAscan-SE -o "$output_folder/tRNAscan-SE_output_$prefix_inside" "$genomic_seq_full_path"
    check_error "tRNAscan-SE for $genomic_seq"
    # echo "Content" > "$output_folder/plik_$prefix_inside"

    echo "Processing for $genomic_seq and $cds_seq completed. It was $counter th pair"
    echo " "
    echo " "
}

# Export variables to parallel function
export -f process_pair
export -f check_error
export genomic_folder_path
export cds_folder_path
export files_per_folder=19
export output_base_folder="output_shorted_pipeline_tAI"
num_pairs=${#file_pairs[@]}

# Calculate number of chunks and folder index management
num_chunks=$(( (num_pairs + files_per_folder - 1) / files_per_folder ))

echo "Processing started"

# Loop to process file pairs
for ((i=0; i<num_pairs; i++)); do
    # Calculate the folder index (1 to 10)
    folder_index=$(( (i / files_per_folder) + 1 ))

    # Ensure the folder exists
    output_folder="output_shorted_pipeline_tAI$folder_index"
    if [ ! -d "$output_folder" ]; then
        mkdir -p "$output_folder"
    fi

    # Process the pair (running process in background)
    process_pair "${file_pairs[$i]}" $((i+1)) "$output_folder" &

    # Wait for the current batch to finish before continuing (every 19 files)
    if (( (i + 1) % files_per_folder == 0 )); then
        wait  # Wait for all background processes to complete
    fi
done

# Ensure all background processes finish
wait

echo "All chunks have been processed."
