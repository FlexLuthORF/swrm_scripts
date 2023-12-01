import csv
import sys
import os
import subprocess
import time

def submit_job(sample_id, data_path, root_path):
    scratch = f"{root_path}/minimap/{sample_id}"
    os.makedirs(scratch, exist_ok=True)
    job_script = "/home/zmvanw01/git_repos/swrm_scripts/zvw/revio_temp.sh"

    # Common directory for SLURM output and error files
    output_dir = os.path.join(root_path, "slurm_outputs")
    os.makedirs(output_dir, exist_ok=True)

    # Construct the SLURM sbatch command
    sbatch_command = [
        "sbatch",
        "--job-name=revio_map",
        f"--output={output_dir}/{sample_id}_%j.out",
        f"--error={output_dir}/{sample_id}_%j.err",
        "--ntasks=1",
        "--cpus-per-task=12",
        job_script,
        data_path,
        scratch
    ]

    # Submit the job
    subprocess.run(sbatch_command)

def get_running_jobs():
    result = subprocess.run(["squeue", "--name=revio_map", "--noheader"], capture_output=True, text=True)
    return len(result.stdout.strip().split('\n'))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <root_path>")
        sys.exit(1)

    root_path = sys.argv[1]
    max_jobs = 10  # Maximum number of concurrent jobs

    with open('/home/watsonlab/project/CW61/STC567/output.tsv', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader)  # Skip header row if present
        for row in reader:
            if len(row) == 2:
                while get_running_jobs() >= max_jobs:
                    time.sleep(60)  # Wait for 60 seconds before checking again
                submit_job(row[0], row[1], root_path)
            else:
                print(f"Invalid row in TSV file: {row}")
