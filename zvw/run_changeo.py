import os
#import glob
from multiprocessing import Pool
import argparse
import subprocess

def changeo(sample_info, outdir):
    sample_id, r1_gz_path = sample_info
    sample_output_dir = os.path.join(outdir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    log_dir = os.path.join(sample_output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    changeo_path= 
    
    r_script_command = [
        "Rscript", 
        "path_to_Rscript", 
        #r1_gz_path
    ]
    
    # Run R script
    subprocess.run(r_script_command, check=True)
    
    # After Rscript has run, execute the CreateGermlines.py command
    # Assuming your R script outputs the necessary file in the correct directory.
    create_germlines_command = [
        "CreateGermlines.py", 
        "-d", os.path.join(sample_output_dir, "changeo", "output_file_from_Rscript.tsv"), 
        "-r", "~/share/germlines/imgt/human/vdj", 
        "-g", "dmask", 
        "--format", "airr", 
        "--cloned"
    ]
    
    # Run CreateGermlines.py
    subprocess.run(create_germlines_command, check=True)



def read_sample_file(file_list_path):
    with open(file_list_path, 'r') as file_list:
        return [line.strip().split() for line in file_list if line.strip()]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some samples.')
    parser.add_argument('--file', required=True, help='Path to the file containing the list of samples.')
    parser.add_argument('--outdir', required=True, help='Desired output directory path.')

    args = parser.parse_args()

    sample_info_list = read_sample_file(args.file)
    with Pool() as pool:
        pool.starmap(process_sample, [(info, args.outdir) for info in sample_info_list])