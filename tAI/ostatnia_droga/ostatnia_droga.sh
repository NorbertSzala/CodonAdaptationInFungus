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
echo "pairs: $file_pairs"   




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
    echo "hej"

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


    # # # # Process CDS sequence using other tools sequentially
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

    # # Clean up intermediate files
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


chunk_size=19
total_pairs=${#file_pairs[@]}

# Liczba chunków
num_chunks=$(( (total_pairs + chunk_size - 1) / chunk_size ))

for ((i=0; i<num_chunks; i++)); do
    # Określenie indeksów dla chunku
    start_index=$((i * chunk_size))
    end_index=$((start_index + chunk_size - 1))

    # Zapewnia, że ostatni chunk nie wyjdzie poza zakres
    if ((end_index >= total_pairs)); then
        end_index=$((total_pairs - 1))
    fi

    # Tworzenie chunku
    chunk=("${file_pairs[@]:$start_index:$((end_index - start_index + 1))}")

    # Zapisz zawartość chunku do pliku
    chunk_file="$output_folder/output_shortened_pipeline_tAI_chunk$((i+1)).txt"
    echo "Zapisuję zawartość chunku $((i+1)) do pliku: $chunk_file"
    
    # Zapisanie zawartości chunku do pliku
    {
        echo "Zawartość chunku $((i+1)):"
        for ((j=0; j<${#chunk[@]}; j++)); do
            echo "${chunk[$j]}"
        done
    } > "$chunk_file" &

    # Uruchomienie operacji na parach w tle w danym chunku
done

# Czekaj na zakończenie wszystkich procesów
wait

echo "Wszystkie operacje zakończone!"
