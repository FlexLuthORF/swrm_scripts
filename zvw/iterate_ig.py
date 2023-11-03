import os
import subprocess
from multiprocessing import Pool
import argparse
from functools import partial


def count_lines(file_path):
    try:
        with open(file_path, "r") as f:
            return sum(1 for _ in f)
    except FileNotFoundError:
        print(f"File {file_path} not found.")
        return 0


def generate_sbatch(num_cores, sample_file, out_dir, job):
    
    if job == 'igk_pipeline':
        taskpernode = '1'
    elif job == 'ig_honda_vcf'or job == 'ig_merge_hifiasm':
        taskpernode = '12' 
    sbatch_content = f"""#!/bin/bash
#SBATCH --job-name={job}
#SBATCH --partition=compute             # Partition (job queue)
#SBATCH --time=01:00:00              # Runtime in D-HH:MM format
#SBATCH --output=job_logs/%j.out     # File to which STDOUT will be written
#SBATCH --error=job_logs/%j.err      # File to which STDERR will be written
#SBATCH --nodes=1-10
#SBATCH --cpus-per-task={taskpernode}
#SBATCH --ntasks-per-node={taskpernode}
#SBATCH --ntasks={num_cores}
echo "Starting job for {sample_file} with {num_cores} cores."

python ./{job}.py --file {sample_file} --outdir {out_dir}

echo "Job finished."

"""

    with open("dynamic_sbatch.sbatch", "w") as f:
        f.write(sbatch_content)
    #print(f"Generated sbatch content:\n{sbatch_content}")
def submit_sbatch():

    result = subprocess.run("sbatch dynamic_sbatch.sbatch", shell=True)
    if result.returncode != 0:
        print(f"Sbatch submission failed with error code {result.returncode}")
    return result.returncode  

def run_IG_for_sample(rpath, sample):
    command = f"bash /path/to/your/run_IG.sh --sample {sample} --outdir {rpath}"
    subprocess.run(command, shell=True)

def run_do_igks(files, rpath):
    with open(files, "r") as file:
        samples = [line.strip().split()[0] for line in file]

    # Run function in parallel
    with Pool() as pool:
        pool.map(partial(run_IG_for_sample, rpath=rpath), samples)

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description='Run IG analysis.')
    argparser.add_argument('--file', '-f', type=str, help='Path to the .txt file containing sample IDs', required=True)
    argparser.add_argument('--outdir', '-o', type=str, help='Output directory', default=os.getcwd())
    
    args = argparser.parse_args()
    rpath = args.outdir
    sample_file_path = args.file
    
    print(sample_file_path)
    run_do_igks(sample_file_path, rpath)
