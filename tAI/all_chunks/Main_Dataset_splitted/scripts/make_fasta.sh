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
    # Dla każdej sekwencji w kolumnach 2, 3 i 4
    for i in 1 2 3; do
        seq_var="seq$i"
        seq_value="${!seq_var}"

        # Jeśli sekwencja istnieje, zapisujemy ją w formacie FASTA
        if [ -n "$seq_value" ]; then
            echo ">$protein_id_seq$i" >> "$output_file"
            echo "$seq_value" >> "$output_file"
        fi
    done
done < "$input_file"
