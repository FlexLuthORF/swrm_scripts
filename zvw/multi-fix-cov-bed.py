import sys
import os
import pandas as pd
from multiprocessing import Pool


def find_element(start, end, ref_bed):
    for _, row in ref_bed.iterrows():
        ref_start = row['start']
        ref_end = row['end']
        if start >= ref_start and end <= ref_end:
            return row['element']
    return "."

def process_sample(sample_id):
    # Read the input BED file for the sample
    input_file = f"{os.getcwd()}/ccs_cov/{sample_id}/{sample_id}_ccs-cov.bed"
    input_bed = pd.read_csv(input_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])

    # Read the reference BED file
    ref_file = "/home/zmvanw01/test-beds/sorted_region.bed"
    ref_bed = pd.read_csv(ref_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'element'])

    # Process each row of the input BED file
    output_bed = input_bed.copy()
    output_bed['element'] = output_bed.apply(lambda row: find_element(row['start'], row['end'], ref_bed), axis=1)
    output_bed['sample_id'] = sample_id

    return output_bed

# Check if the fofn file is provided as an argument
if len(sys.argv) < 2:
    print("Please provide the fofn file as an argument.")
    sys.exit(1)

# Read the fofn file from the command-line argument
fofn_file = sys.argv[1]
with open(fofn_file, "r") as file:
    sample_ids = [line.strip().split("\t")[0] for line in file]

# Create a multiprocessing pool with 12 threads
pool = Pool(processes=12)

# Process each sample in parallel
results = pool.map(process_sample, sample_ids)

pool.close()
pool.join()

# Concatenate all the processed DataFrames into a single DataFrame
merged_bed = pd.concat(results, ignore_index=True)

# Write the merged DataFrame to a master TSV file
output_file = f"{os.getcwd()}/ccs_cov/master_output.tsv"
merged_bed.to_csv(output_file, sep='\t', header=True, index=False)

print(f"Merged output written to {output_file}")