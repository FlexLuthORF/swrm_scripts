import argparse
import os
import subprocess
from multiprocessing import Pool

# Function to process each sample
def process_sample(sample_id, root_dir):
    # Check if the TSV file exists
    tsv_path = os.path.join(root_dir, "presto", sample_id, "S5_db-pass_clone-pass_germ-pass.tsv")
    #print(tsv_path)
    if os.path.exists(tsv_path):
        # Call the R script with sample_id and root_dir as arguments
        subprocess.run(["Rscript", "/home/zmvanw01/git_repos/swrm_scripts/zvw/counts.R", sample_id, root_dir])

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process samples from a text file")
    parser.add_argument("--dir", required=True, help="Root directory")
    parser.add_argument("--file", required=True, help="Path to the samples text file")
    args = parser.parse_args()

    # Read sample IDs from the text file
    with open(args.file, "r") as file:
        sample_ids = [line.strip().split()[0] for line in file]

    # Create a pool of workers for parallel processing
    with Pool() as pool:
        # Use starmap to pass sample_id and root_dir to the process_sample function
        pool.starmap(process_sample, [(sample_id, args.dir) for sample_id in sample_ids])
