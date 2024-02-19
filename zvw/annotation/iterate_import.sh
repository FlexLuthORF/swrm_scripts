#!/bin/bash

# Check if an argument was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 path/to/your/fofn.tsv"
    exit 1
fi

# Use the first argument as the path to your fofn.tsv file
fofn_path="$1"

# Directory where the CSV files are located
csv_dir=$(pwd)

while IFS=$'\t' read -r sampleID _; do
    python /home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/import_from_assemblies.py IGK chr2 "+-" "${csv_dir}/${sampleID}_genes.csv" /home/zmvanw01/git_repos/VDJbase-Genomics/results/IGK/genome/chr2.fasta /home/zmvanw01/git_repos/VDJbase-Genomics/results/IGK/coords ~/git_repos/VDJbase-Genomics/ref/ "${csv_dir}/${sampleID}_genes_imported.csv"
done < "$fofn_path"
