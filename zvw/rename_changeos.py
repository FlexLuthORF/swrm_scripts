import pandas as pd
import sys
import os
import numpy as np

# Check if command line argument is provided
if len(sys.argv) < 2:
    print("Usage: python rename_tsv.py <path_to_sampleID.txt>")
    sys.exit(1)

sample_ids_file = sys.argv[1]
guide_file = "/home/zmvanw01/projects/repseq-files/collapse_files/merged_rename_corrected.tsv"
missing_samples_file = "missing_samples.txt"

# Special case sample IDs
special_sample_ids = {
    "211170201C", "211281405C", "220173102C", "2202414003", "2202415007", "2202423018", "2202424002", 
    "220271402C", "220272803C", "220370304C", "220372201C", "220372302C", "220373103C", "220380302C", 
    "220380901C", "220381703C", "2204419003", "220470501C", "220472101C", "220481803C", "220482104C", 
    "220770702C", "2203414011", "2208409002", "220881602C", "220881704C", "2209401001", "2209420003", 
    "220980104C", "220981402C", "220981504C", "2210410005", "2210411018", "2211414001", "190980901C", 
    "190981003C", "190981502C", "191080609C", "2203414003", "220372304C", "220372802C", "220373102C", 
    "220382303C", "220383004C", "220471901C", "220472801C", "220580504C", "220581202C", "220670601C", 
    "220672801C", "220680803C", "220680902C", "220681401C", "220682301C", "220770701C", "220771101C", 
    "220772601C", "220781702C", "220782501C", "220782502C", "220782701C", "2208411000", "220881003C", 
    "220882403C", "2209430000", "220971201C", "220971301C", "220972003C", "220972101C", "220980801C", 
    "220982201C", "220982202C", "2210428001", "221071801C", "2211429004", "221170702C", "221170803C", 
    "221170901C", "221171504C", "2301503003", "230170402C", "230170502C", "230170902C", "2212405008", 
    "2301423007", "2301425004", "2301826101", "2301826110", "2209413009", "2210419015", "230271302C", 
    "2302807101", "2302809105", "2303822105", "2303824104"
}

# Read the guide file
guide_df = pd.read_csv(guide_file, sep='\t')

# Define a function to rename genes
def rename_gene(value, genes, change_to):
    if isinstance(value, str) and not pd.isna(value):
        if any(gene in value for gene in genes):
            return change_to
    return value

# Process each sample ID
with open(sample_ids_file, 'r') as file, open(missing_samples_file, 'w') as missing_file:
    for line in file:
        sample_id = line.split()[0].strip()
        tsv_file = f"./presto/{sample_id}/changeo/{sample_id}_merged-changeo.tsv"
        out_file = f"./presto/{sample_id}/changeo/{sample_id}_merged-changeo-rename.tsv"

        if not os.path.exists(tsv_file):
            missing_file.write(sample_id + '\n')
            continue

        try:
            df = pd.read_csv(tsv_file, sep='\t', low_memory=False)

            for _, row in guide_df.iterrows():
                # Check if 'genes' value is not NaN
                if pd.notna(row['genes']):
                    genes = row['genes'].split('_')
                    change_to = row['change_to']
                    
                    # Skip special case if the sample ID matches
                    if genes == ["IGKV1-13", "IGKV1D-13"] and sample_id in special_sample_ids:
                        continue

                    for gene in genes:
                        column = f"{gene[3].lower()}_call"
                        if column in df.columns:
                            df[column] = df[column].apply(rename_gene, args=(genes, change_to))
                else:
                    # Skip and report NaN genes value
                    print(f"Skipping NaN genes value in row: {row['genes']}")

            df.to_csv(out_file, sep='\t', index=False)
        except FileNotFoundError:
            print(f"File not found: {tsv_file}")
        except Exception as e:
            print(f"Error processing {sample_id}: {e}")

print("Processing complete.")
