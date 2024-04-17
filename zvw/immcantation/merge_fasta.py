
import os
import glob
import shutil

input_dir = "."  # Replace with your actual directory
fasta_filename = "S5-final_total_collapse-unique_atleast-2_reheader.fasta"
#output_dir, f"{sample_id}.fasta"), 'w'

for subdir in glob.glob(os.path.join(input_dir, "*_*")):
    sample_id = subdir.split("_")[0]
    output_dir = os.path.join(input_dir, f"{sample_id}_merged")
    os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

    fasta_files = glob.glob(os.path.join(subdir, fasta_filename))  # Find specific file
    output_file = os.path.join(output_dir, fasta_filename)

    with open(output_file, 'w') as outfile:
        for fasta_file in fasta_files:
            with open(fasta_file, 'r') as infile:
                shutil.copyfileobj(infile, outfile) 
