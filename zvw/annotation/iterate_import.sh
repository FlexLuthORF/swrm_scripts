#!/bin/bash

# Check if an argument was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 path/to/your/fofn.tsv"
    exit 1
fi

# Use the first argument as the path to your fofn.tsv file
fofn_path="$1"

# Base directory where the CSV files are located
base_dir=$(pwd)

# Array of directories to loop through
dirs=("igh" "igk" "igl")

# Loop through each directory
for dir in "${dirs[@]}"; do
    # Set the command options based on the directory
    case $dir in
        igh)
            options=("IGH" "igh" '\"+-\"')
            ;;
        igk)
            options=("IGK" "chr2" '\"+-\"')
            ;;
        igl)
            options=("IGL" "chr22" "+")
            ;;
    esac

    # Directory for current iteration
    csv_dir="${base_dir}/${dir}"

    # Loop through the fofn.tsv file
    while IFS=$'\t' read -r sampleID _; do
        # Use eval to correctly handle the inclusion of quotes in options
        eval python /home/zmvanw01/git_repos/swrm_scripts/zvw/annotation/import_from_assemblies.py "${options[0]}" "${options[1]}" "${options[2]}" "${csv_dir}/${sampleID}_genes.csv" /home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta /home/zmvanw01/git_repos/immune_receptor_genomics/current /home/zmvanw01/git_repos/VDJbase-Genomics/ref/ "${csv_dir}/${sampleID}_genes_imported.csv"
    done < "$fofn_path"
done
