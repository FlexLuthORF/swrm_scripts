import csv
import sys
import os
from multiprocessing import Pool
import subprocess

def run_script(sample_id, data_path):
    scratch = f"{root_path}/minimap/{sample_id}"

    # Check if scratch exists, if not, create it
    if not os.path.exists(scratch):
        os.makedirs(scratch)

    subprocess.run(["/home/zmvanw01/git_repos/swrm_scripts/zvw/revio_temp.sh", data_path, scratch], check=True)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <root_path>")
        sys.exit(1)

    root_path = sys.argv[1]
    samples = []

    with open('/home/watsonlab/project/CW61/STC567/output.tsv', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header row if present
        for row in reader:
            if len(row) == 2:
                samples.append(tuple(row))
            else:
                print(f"Invalid row in TSV file: {row}")

    with Pool() as pool:
        pool.starmap(run_script, samples)
