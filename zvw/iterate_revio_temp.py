import csv
import sys
import os
import subprocess
import time
import glob

def submit_job(sample_id, data_path, root_path):
    scratch = f"{root_path}/minimap/{sample_id}"
    os.makedirs(scratch, exist_ok=True)
    job_script = "/home/zmvanw01/git_repos/swrm_scripts/zvw/revio_temp.sh"

    output_dir = os.path.join(root_path, "slurm_outputs")
    os.makedirs(output_dir, exist_ok=True)

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

    subprocess.run(sbatch_command)

def get_running_jobs():
    result = subprocess.run(["squeue", "--name=revio_map", "--noheader"], capture_output=True, text=True)
    return len(result.stdout.strip().split('\n'))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <root_path> <csv_file>")
        sys.exit(1)

    root_path = sys.argv[1]
    csv_file = sys.argv[2]
    max_jobs = 10

    with open(csv_file, newline='') as file:
        reader = csv.reader(file, delimiter=',')
        next(reader)
        for row in reader:
            if len(row) > 0:
                sample_id = row[0]
                data_path_pattern = f"{root_path}/{sample_id}/{sample_id}*.bam"
                data_paths = glob.glob(data_path_pattern)
                if data_paths:
                    data_path = data_paths[0]
                    while get_running_jobs() >= max_jobs:
                        time.sleep(60)
                    submit_job(sample_id, data_path, root_path)
                else:
                    print(f"No data found for {sample_id}")
            else:
                print("Invalid row in CSV file")
