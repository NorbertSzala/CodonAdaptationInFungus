#!/bin/bash

# Sprawdzamy, czy zostały podane odpowiednie argumenty
if [ $# -ne 2 ]; then
    echo "Użycie: $0 <plik_wejściowy.tsv> <plik_wyjściowy.fasta>"
    exit 1
fi

# Przypisujemy argumenty do zmiennych
input_file="$1"
output_file="$2"

# Sprawdzamy, czy plik wejściowy istnieje
if [ ! -f "$input_file" ]; then
    echo "Błąd: Plik wejściowy $input_file nie istnieje."
    exit 1
fi

# Przechodzimy przez każdy wiersz w pliku wejściowym
while IFS=$'\t' read -r protein_id seq1 seq2 seq3; do
    # Jeśli istnieje sekwencja 1, zapisujemy ją w formacie FASTA
    if [ -n "$seq1" ]; then
        echo ">$protein_id_seq1" >> "$output_file"
        echo "$seq1" >> "$output_file"
    fi
    # Jeśli istnieje sekwencja 2, zapisujemy ją w formacie FASTA
    if [ -n "$seq2" ]; then
        echo ">$protein_id_seq2" >> "$output_file"
        echo "$seq2" >> "$output_file"
    fi
    if [ -n "$seq3" ]; then
        echo ">$protein_id_seq3" >> "$output_file"
        echo "$seq3" >> "$output_file"
    fi
done < "$input_file"

