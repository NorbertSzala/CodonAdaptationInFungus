#!/bin/bash

# Zmienna przechowująca wynikowy plik
output_file1="concatenated_each_protein_file.txt"
output_file2="concatenated_summary_protein_file.txt"
output_file3="concatenated_each_protein_file_shorted.txt"
output_file4="concatenated_summary_protein_file_shorted.txt"
# Zainicjalizuj plik wynikowy (jeśli już istnieje, zostanie nadpisany)
> "$output_file1"
> "$output_file2"


# Iteracja przez foldery chunk1, chunk2, ..., chunk10
for i in {1..10}; do
  # Ścieżka do folderu
  folder="/home/norbert_s/tAI/all_chunks/chunk$i/tAI_folder"
  echo "$folder"
  
#   # Iteracja przez pliki w folderze
  for file in "$folder"/each_protein_tAI_*.txt; do
    # Pobierz nazwę pliku bez ścieżki
    filename=$(basename "$file")
    
    # Wyciągnij prefix (13 ostatnich znaków przed rozszerzeniem .txt)
    prefix=$(echo "$filename" | sed -E 's/.*(GCA_[0-9]{9})\.txt/\1/')

    # Dodaj prefix jako pierwszą kolumnę i połącz plik z plikiem wynikowym
    awk -v prefix="$prefix" '{print prefix"\t"$0}' "$file" >> "$output_file1"
    
  done
  sed -i '1s/^[^ ]*/assembly_ID\t"title"\t"Nc"\t"GC3s"\t"GC"\t"L_aa"\t"AA_len_bin"\t"tAI"/' "$output_file1"
  
done
awk '{print $2, $NF}' "$output_file1" > "$output_file3"


# Iteracja przez foldery chunk1, chunk2, ..., chunk10
for i in {1..10}; do
  # Ścieżka do folderu
  folder="/home/norbert_s/tAI/all_chunks/chunk$i/tAI_folder"
  echo "$folder"
  
  # Iteracja przez pliki w folderze
  for file in "$folder"/summary_tAI_*.txt; do
    # Pobierz nazwę pliku bez ścieżki
    filename=$(basename "$file")
    
    # Wyciągnij prefix (13 ostatnich znaków przed rozszerzeniem .txt)
    prefix=$(echo "$filename" | sed -E 's/.*(GCA_[0-9]{9})\.txt/\1/')

    echo "$prefix"
    # Pomijaj pierwszy wiersz i dodaj prefix jako pierwszą kolumnę
    awk -v prefix="$prefix" 'NR > 1 {print prefix"\t"$0}' "$file" >> "$output_file2"

  done
  sed -i '1s/^[^ ]*/assembly_ID\tpath\t"tAI_mean"\t"GC_mean"\t"L_aa_mean"\t"pearson"\t"correlation_tAI_NC_GC3"/' "$output_file2"

done
awk '{print $1, $3}' "$output_file2" > "$output_file4"