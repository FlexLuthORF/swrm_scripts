import sys
import pandas as pd

def find_element(start, end, ref_bed):
    for _, row in ref_bed.iterrows():
        ref_start = row['start']
        ref_end = row['end']
        if start >= ref_start and end <= ref_end:
            return row['element']
    return "."

# Check if the input BED file is provided as an argument
if len(sys.argv) < 2:
    print("Please provide the input BED file as an argument.")
    sys.exit(1)

# Read the input BED file from the command-line argument
input_file = sys.argv[1]
input_bed = pd.read_csv(input_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])

# Read the reference BED file
ref_file = "/home/zmvanw01/test-beds/sorted_region.bed"
ref_bed = pd.read_csv(ref_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'element'])

# Process each row of the input BED file
output_bed = input_bed.copy()
output_bed['element'] = output_bed.apply(lambda row: find_element(row['start'], row['end'], ref_bed), axis=1)

# Write the output BED file
output_file = "output.bed"
output_bed.to_csv(output_file, sep='\t', header=False, index=False)

print("Output written to", output_file)