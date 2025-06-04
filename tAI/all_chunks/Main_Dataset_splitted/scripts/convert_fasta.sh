#!/bin/bash

if [ -z "$1" ]; then
    echo "Proszę podać plik wejściowy jako argument."
    exit 1
fi

input_file="$1"

if [ ! -f "$input_file" ]; then
    echo "Plik '$input_file' nie istnieje."
    exit 1
fi


if [ -z "$2" ]; then
    echo "Proszę podać plik wyjściowy jako drugi argument."
    exit 1
fi

output_file="$2"


while IFS=$'\t' read -r id seq1 seq2 seq3; do

    if [[ -n "$id" ]]; then
        # Numeracja sekwencji
        seq_number=1
	echo "$id"

        if [[ -n "$seq1" ]]; then
            echo ">$id_seq$seq_number" >> "$output_file"
            echo "$seq1" >> "$output_file"
            ((seq_number++))
        fi
        if [[ -n "$seq2" ]]; then
            echo ">$id_seq$seq_number" >> "$output_file"
            echo "$seq2" >> "$output_file"
            ((seq_number++))
        fi
        if [[ -n "$seq3" ]]; then
            echo ">$id_seq$seq_number" >> "$output_file"
            echo "$seq3" >> "$output_file"
            ((seq_number++))
        fi
    fi
done < "$input_file"


