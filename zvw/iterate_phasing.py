import os
import argparse
import subprocess

def count_subdirectories(root_dir):
    return len([name for name in os.listdir(root_dir) if os.path.isdir(os.path.join(root_dir, name))])

def create_sbatch_file(root_dir, script_path, sample_count):
    sbatch_content = f"""#!/bin/bash
#SBATCH --job-name=phase
#SBATCH --partition=compute             # Partition (job queue)
#SBATCH --output=job_logs/%j.out     # File to which STDOUT will be written
#SBATCH --error=job_logs/%j.err      # File to which STDERR will be written
#SBATCH --nodes=10
#SBATCH --cpus-per-task=12
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks={sample_count}
echo "Starting job for phasing with {sample_count} samples."
python {script_path}
"""

    sbatch_file = os.path.join(root_dir, "submit_job.sbatch")
    with open(sbatch_file, "w") as file:
        file.write(sbatch_content)
    
    return sbatch_file

def main():
    parser = argparse.ArgumentParser(description='Generate and submit SBATCH file for sample processing.')
    parser.add_argument('root_dir', help='Root directory containing subdirectories with .bam files')
    args = parser.parse_args()

    sample_count = count_subdirectories(args.root_dir)
    script_path = '/home/zmvanw01/git_repos/swrm_scripts/zvw/run_phasing.py'  # Replace with actual path of the first script
    sbatch_file = create_sbatch_file(args.root_dir, script_path, sample_count)

    # Submit the SBATCH job without waiting for completion
    subprocess.Popen(["sbatch", sbatch_file])

if __name__ == "__main__":
    main()
