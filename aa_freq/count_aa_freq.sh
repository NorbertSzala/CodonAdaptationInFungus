#!/bin/bash

#check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
	echo "Usage $0 [input file] [output file]"
	exit 1
fi

#Assign command-line arguments to variables
input_file=$1
output_file=$2

#check if the input file exist
if [ ! -f "$input_file" ]; then
	echo "Error: Input file '$input_file' does not exist."
	exit 1
fi

#iterate through all aminoacids and write its count to file
amino_acids=("A Ala" "C Cys" "D Asp" "E Glu" "F Phe" "G Gly" "H His" "I Ile" "K Lys" "L Leu" "M Met" "N Asn" "P Pro" "Q Gln" "R Arg" "S Ser" "T Thr" "V Val" "W Trp" "Y Tyr")

sum=0

for aa in "${amino_acids[@]}"; do
	one_letter=$(echo $aa | cut -d ' ' -f 1) # cut -d ' ' sets a space as delimiter/separator. -f 1 chooses first newly created column
	three_letter=$(echo $aa | cut -d ' ' -f 2)
	count=$(grep -v ">" "$1" | grep -o "$one_letter" | wc -l)
	sum=$((sum+count))
#	echo "$one_letter $three_letter: $count" >> "$2"
	printf "%s\t%s\t%d\n" "$three_letter" "$one_letter" "$count" >> "$2" #%s means string \t tabulator %d number
done

echo "Total count: $sum"


