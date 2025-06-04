#!/bin/bash

# ================================================
# Script to Extract and Merge tAI Data from R Files
#
# This bash script processes tAI-related data by extracting
# and merging results from multiple R-generated output files.
# It creates concatenated tables with both detailed and summarized
# information for further analysis.
#
# Key Functionality:
# 1. **Output Files**:
#    - `output_file1`: Concatenated protein-level data with tAI and related values.
#    - `output_file2`: Concatenated summary-level data with tAI means and other statistics.
#    - `output_file3`: Shortened version of `output_file1` with essential values.
#    - `output_file4`: Shortened version of `output_file2` with key statistics.
#
# 2. **Data Extraction and Merging**:
#    - The script iterates through directories `chunk1` to `chunk10`, processing files in `tAI_folder`.
#    - It extracts the assembly ID from filenames, adds it as a new column, and appends the data to the output files.
#
# 3. **File Formatting**:
#    - Adds headers to the output files.
#    - Creates shortened versions with only essential columns.
#
# 4. **Final Output**:
#    - Generates four output files: `output_file1`, `output_file2`, `output_file3`, and `output_file4`.
#
# Dependencies:
# - Assumes that input R files with tAI data are already generated.
#
# Usage:
# To run the script, execute the following command:
# bash data_merging_script.sh
# ================================================


# Variable to store the output file
output_file1="concatenated_each_protein_file.txt"
output_file2="concatenated_summary_protein_file.txt"
output_file3="concatenated_each_protein_file_shorted.txt"
output_file4="concatenated_summary_protein_file_shorted.txt"


# Iterating through chunk1, chunk2, ..., chunk10 folders
for i in {1..10}; do
  # Path to the folderu
  folder="/home/norbert_s/tAI/all_chunks/chunk$i/tAI_folder"
  echo "$folder"
  
# Iterating through files in the folder
  for file in "$folder"/each_protein_tAI_*.txt; do
     # Get the filename without the path
    filename=$(basename "$file")
    
    # Extract the prefix (15 last characters before the .txt extension)
    prefix=$(echo "$filename" | sed -E 's/.*(GCA_[0-9]{9})\.txt/\1/')

    # Add prefix as the first column and merge the file with the result file
    awk -v prefix="$prefix" '{print prefix"\t"$0}' "$file" >> "$output_file1"
  done

  sed -i '1s/^[^ ]*/assembly_ID\tprotein_id\tNc\tGC3s\tGC\tL_aa\tAA_len_bin\ttAI/' "$output_file1"
  
done

awk '{print $1 "\t" $2 "\t" $NF}' "$output_file1" > "$output_file3"


# Iterating through chunk1, chunk2, ..., chunk10 folders
for i in {1..10}; do
  # Path to the folder
  folder="/home/norbert_s/tAI/all_chunks/chunk$i/tAI_folder"
  echo "$folder"

  # Iterating through files in the folder
  for file in "$folder"/summary_tAI_*.txt; do
    # Get the filename without the path
    filename=$(basename "$file")

	cut -f2- "$file" > temp_file


    # Extract the prefix (15 last characters before the .txt extension)
    prefix=$(echo "$filename" | sed -E 's/.*(GCA_[0-9]{9})\.txt/\1/')

    awk -v prefix="$prefix" 'NR>1 {print prefix"\t"$0}' temp_file >> "$output_file2"


  done
done
sed -i "1i $(echo -e 'assembly_ID\ttAI_mean\tGC_mean\tL_aa_mean\tpearson\tcorrelation_tAI_NC_GC3')" "$output_file2"

awk '{print $1 "\t" $2}' "$output_file2" > "$output_file4"

echo "Done."