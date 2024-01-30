import os
import subprocess
import time
import argparse

def count_running_jobs():
    """Count the number of running sbatch jobs."""
    result = subprocess.run(["squeue", "-u", os.environ["USER"]], capture_output=True, text=True)
    return result.stdout.count("\n") - 1  # Subtracting header line

def submit_job(sample_id, r1_gz_path, script_path):
    """Submit a job using sbatch with inline command and custom log files."""
    log_dir = "/home/zmvanw01/projects/CW50/CW50_NS_pool3"
    stdout_log = os.path.join(log_dir, f"{sample_id}_stdout.log")
    stderr_log = os.path.join(log_dir, f"{sample_id}_stderr.log")

    command = (
        f"sbatch "
        f"--output='{stdout_log}' "
        f"--error='{stderr_log}' "
        f"--export=SAMPLE_ID='{sample_id}',R1_GZ_PATH='{r1_gz_path}' "
        f"--wrap='python {script_path}'"
    )
    subprocess.run(command, shell=True)

def read_sample_file(file_list_path):
    """Read sample information from the file."""
    with open(file_list_path, 'r') as file_list:
        return [line.strip().split() for line in file_list if line.strip()]

def main(file_list_path):
    processing_script_path = '/home/zmvanw01/git_repos/swrm_scripts/zvw/immcantation/run_presto.py'
    max_jobs = 15

    sample_info_list = read_sample_file(file_list_path)

    for sample_id, r1_gz_path in sample_info_list:
        while count_running_jobs() >= max_jobs:
            time.sleep(60)  # Wait for 60 seconds before checking again
        submit_job(sample_id, r1_gz_path, processing_script_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process file paths.')
    parser.add_argument('file_list_path', type=str, help='Path to the sample file list')
    args = parser.parse_args()

    main(args.file_list_path)
