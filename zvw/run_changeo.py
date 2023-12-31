import os
#import glob
from multiprocessing import Pool
import argparse
import subprocess

def changeo(sample_info, outdir):
    sample_id, r1_gz_path = sample_info
    sample_output_dir = os.path.join(outdir, sample_id)
    sample_changeo_db = os.path.join(sample_output_dir,"S5_db-pass.tsv")
    os.makedirs(sample_output_dir, exist_ok=True)
    log_dir = os.path.join(sample_output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)
    print(sample_changeo_db)

    r_script_command = [
        "Rscript", 
        "/home/zmvanw01/git_repos/swrm_scripts/zvw/defineclones.R", 
        sample_changeo_db
    ]
    
    if not os.path.exists(os.path.join(sample_output_dir, "S5_db-pass_clone-pass.tsv")):
        subprocess.run(r_script_command, check=True)
    
    
    create_germlines_command = [
        "/home/zmvanw01/anaconda3/envs/immcantation/bin/python", "/home/zmvanw01/.local/bin/CreateGermlines.py", 
        "-d", os.path.join(sample_output_dir, "S5_db-pass_clone-pass.tsv"), 
        "-r", "/home/zmvanw01/share/germlines/imgt/human/vdj", 
        "-g", "dmask", 
        "--format", "airr", 
        "--cloned"
    ]
    
    if not os.path.exists(os.path.join(sample_output_dir, "S5_db-pass_clone-pass-germ_pass.tsv")):
        print(create_germlines_command)
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
        pool.starmap(changeo, [(info, args.outdir) for info in sample_info_list])
