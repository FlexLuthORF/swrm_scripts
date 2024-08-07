import subprocess
import argparse
import csv
import xml.etree.ElementTree as ET
import os
import shutil

# Argument parser setup
parser = argparse.ArgumentParser(description='Run SCP commands and extract information.')
parser.add_argument('base_directory', help='Base directory path for SCP command.')
args = parser.parse_args()

# SSH Key path
ssh_key_path = "/home/zmvanw01/.ssh/id_rsa_hydra"

# SSH command to list directories
ssh_command = f"ssh -i {ssh_key_path} zachvanwinkle@136.165.158.10 'ls -d {args.base_directory}/*/'"
result = subprocess.run(ssh_command, shell=True, capture_output=True, text=True)

# Parse output to get directory names
directories = result.stdout.strip().split("\n")

# Initialize TSV writer
tsv_file = open('info.tsv', 'w', newline='')
tsv_writer = csv.writer(tsv_file, delimiter='\t')

# Special case: handle files in base_directory
scp_base = [
    f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{args.base_directory}/*.bam .",
    f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{args.base_directory}/*.pbi .",
    f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{args.base_directory}/*.xml ."
]
for cmd in scp_base:
    subprocess.run(cmd, shell=True)

# Function to rename and move files
def rename_and_move_files(base_filename, bio_sample_name):
    file_types = ['bam', 'bam.pbi', 'consensusreadset.xml']
    paths = {}
    for f_type in file_types:
        original_file = f"{base_filename}.{f_type}"
        new_file = f"{bio_sample_name}_{base_filename}.{f_type}"
        if os.path.exists(original_file):
            os.rename(original_file, new_file)
            new_dir = os.path.abspath(bio_sample_name)
            if not os.path.exists(new_dir):
                os.makedirs(new_dir)
            new_path = os.path.join(new_dir, new_file)
            shutil.move(new_file, new_path)
            if f_type == 'bam':
                paths['bam'] = new_path
    return paths

# Iterate over directories
for dir_path in directories:
    # SCP commands
    scp_bam = f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{dir_path}/*bam ."
    scp_bam_pbi = f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{dir_path}/*bam.pbi ."
    scp_xml = f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{dir_path}/*.xml ."
    
    # Execute SCP commands
    subprocess.run(scp_bam, shell=True)
    subprocess.run(scp_bam_pbi, shell=True)
    subprocess.run(scp_xml, shell=True)

    # Parse XML files and rename & move associated files
    for xml_file in os.listdir('.'):
        if xml_file.endswith('.xml'):
            base_filename = xml_file.rsplit('.consensusreadset.xml', 1)[0]
            try:
                tree = ET.parse(xml_file)
                root = tree.getroot()
                bio_sample_name = root.find(".//{http://pacificbiosciences.com/PacBioSampleInfo.xsd}BioSample").get('Name')
                paths = rename_and_move_files(base_filename, bio_sample_name)
                if 'bam' in paths:
                    tsv_writer.writerow([bio_sample_name, paths['bam']])
            except Exception as e:
                bio_sample_name = "N/A"

# Close TSV file
tsv_file.close()
