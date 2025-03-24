#!/bin/bash

#this script changes fasta headers on needs CodonW. Codonw uses as unique name only first 20 characters from fasta header



input="/db/183_fungi_mai_2022/cds"
output="/home/norbert_s/tAI/all_chunks/cds"

if [ ! -d "$output" ]; then
  mkdir -p "$output"
fi

echo "baza 1"
if [ -d "$input" ]; then
  for file in "$input"/*; do
    if [ -f "$file" ]; then
      cp "$file" "$output/$(basename "$file")"
    fi
  done
fi



for input_file in "$output"/*.fna; do
  if [ -f "$input_file" ]; then
    # Używamy sed z odpowiednim wyrażeniem regularnym
#	protein_id=$(grep "^>" "$input_file" | sed -n 's/.*\[protein_id=\([^]]*\)\].*/\1/p')

    sed -i -E 's/^>.*\[protein_id=([A-Za-z0-9.]+)\].*/>\1/' "$input_file"
#    sed -i -E 's/^[[:space:]]+|[[:space:]]+$//g' "$input_file"


    echo "Plik $input_file przetworzony i nadpisany"
  fi
done