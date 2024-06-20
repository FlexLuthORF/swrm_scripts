import subprocess
import argparse
import csv
import xml.etree.ElementTree as ET
import os
import shutil

# Argument parser setup
parser = argparse.ArgumentParser(description='Run SCP commands and extract information.')
parser.add_argument('base_directory', help='Base directory path for SCP command.')
parser.add_argument('--failed', action='store_true', help='Include failed files in the download.')
args = parser.parse_args()

# SSH Key path
ssh_key_path = "/home/zmvanw01/.ssh/id_rsa_hydra"

# Function to check and download files based on the 'fail' keyword
def check_and_download_files(dir_path, file_types):
    for file_type in file_types:
        # List files on the remote server
        list_command = f"ssh -i {ssh_key_path} zachvanwinkle@136.165.158.135 'find {dir_path} -type f -name \"*.{file_type}\"'"
        result = subprocess.run(list_command, shell=True, capture_output=True, text=True)
        files = result.stdout.strip().split("\n")
        
        # Filter and download files
        for file in files:
            if 'fail' in file.lower() and not args.failed:
                continue
            download_command = f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.135:{file} ."
            subprocess.run(download_command, shell=True)

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
ssh_command = f"ssh -i {ssh_key_path} zachvanwinkle@136.165.158.135 'ls -d {args.base_directory}/*/'"
result = subprocess.run(ssh_command, shell=True, capture_output=True, text=True)
directories = result.stdout.strip().split("\n")

# Initialize TSV writer
tsv_file = open('info.tsv', 'w', newline='')
tsv_writer = csv.writer(tsv_file, delimiter='\t')

for dir_path in directories:
    check_and_download_files(dir_path, ['bam', 'bam.pbi', 'xml'])

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
                print(f"Error processing {xml_file}: {e}")

# Close TSV file
tsv_file.close()
os.system("sed -i 's/\\r//g' info.tsv")
