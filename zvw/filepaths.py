import os
import sys

# Check if the script is given the required argument
if len(sys.argv) < 2:
    print("Usage: script.py <path_to_samples_file>")
    sys.exit(1)

# Path to the samples.txt file (taken from the first command-line argument)
samples_file = sys.argv[1]

# Path to the base directory where the sample folders are located
# Assuming the base directory is the parent directory of samples_file
base_dir = os.path.dirname(samples_file)

# Output file path
output_file = "samples_with_paths.txt"

# Function to find the R1 fastq.gz file for a given sample ID
def find_r1_file(sample_id):
    sample_dir = os.path.join(base_dir, sample_id)
    for root, dirs, files in os.walk(sample_dir):
        for file in files:
            if file.endswith("R1_001.fastq.gz"):
                return os.path.join(root, file)
    return "File not found"

# Reading sample IDs and finding their corresponding R1 file paths
with open(samples_file, 'r') as samples, open(output_file, 'w') as output:
    for sample_id in samples:
        sample_id = sample_id.strip()
        r1_file_path = find_r1_file(sample_id)
        output.write(f"{sample_id}\t{r1_file_path}\n")

print("Process completed. Check the file:", output_file)
