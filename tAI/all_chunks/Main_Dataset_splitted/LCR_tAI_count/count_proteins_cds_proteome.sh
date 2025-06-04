#!/bin/bash


proteome_path='/db/183_fungi_mai_2022/proteomes'
cds_path='/db/183_fungi_mai_2022/cds'

#make sure paths exist
if [ ! -d "$proteome_path" ] || [ ! -d "$cds_path" ]; then
    echo "One or both required paths do not exist:" >&2
    echo "  Proteome path: $proteome_path" >&2
    echo "  CDS path: $cds_path" >&2
    exit 1
fi


for proteom in "$proteome_path"/GCA_*.faa; do
    file_id=$(basename "$proteom" _protein.faa)

    cds_file="$cds_path/${file_id}cds_from_genomic.fna"

    if [ -f "$cds_file" ]; then
        echo "Processing pair: ${file_id}"
        # your processing command here
        













    else
        echo "Missing CDS file for $file_id" >&2
    fi
done
