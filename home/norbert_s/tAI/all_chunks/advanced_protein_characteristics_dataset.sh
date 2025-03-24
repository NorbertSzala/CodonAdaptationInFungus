#!/bin/bash

# Plik wejściowy i wyjściowy
input_file="shorted_main_dataset.tsv"  # Zamień na swój plik wejściowy
output_file="temp_result.txt"  # Zamień na swój plik wyjściowy


# Funkcja do obliczania binu tAI (0-19)
get_tai_bin() {
    echo $(echo "$1" | awk '{print int($1 * 20)}')
}

# Tworzymy plik wyjściowy
> "$output_file"

# Przetwarzamy plik wejściowy
while read -r line; do
    # Rozdzielamy dane na zmienne (kolumny 3, 4, 14)
    protein_length=$(echo "$line" | cut -d' ' -f3)
    domain=$(echo "$line" | cut -d' ' -f4)
    tAI=$(echo "$line" | cut -d' ' -f14)

    # Obliczanie binu tAI
    tai_bin=$(get_tai_bin "$tAI")

    # Dodajemy białko do odpowiednich zbiorów w zależności od domeny
    if [ "$domain" -eq 1 ]; then
        echo "$tai_bin $protein_length 1" >> "$output_file"
    else
        echo "$tai_bin $protein_length 0" >> "$output_file"
    fi
done < "$input_file"  # Wczytujemy dane z pliku wejściowego

# Obliczamy średnie długości białek w każdym binie
awk '
{
    domain=$3; tai_bin=$1; length=$2;
    count[tai_bin, domain]++; sum[tai_bin, domain]+=length;
}
END {
    for (key in count) {
        split(key, idx, SUBSEP);
        tai_bin=idx[1]; domain=idx[2];
        print domain, tai_bin, sum[tai_bin, domain] / count[tai_bin, domain];
    }
}
' "$output_file" > "$output_file"
