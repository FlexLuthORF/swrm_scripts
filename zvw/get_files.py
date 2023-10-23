import subprocess
import argparse
import csv
import xml.etree.ElementTree as ET

# Argument parser setup
parser = argparse.ArgumentParser(description='Run SCP commands and extract information.')
parser.add_argument('base_directory', help='Base directory path for SCP command.')
args = parser.parse_args()

# SSH Key path
ssh_key_path = "~/.ssh/id_rsa_laptop_hydra"

# SSH command to list directories
ssh_command = f"ssh -i {ssh_key_path} zachvanwinkle@136.165.158.10 'ls -d {args.base_directory}/*/'"
result = subprocess.run(ssh_command, shell=True, capture_output=True, text=True)

# Parse output to get directory names
directories = result.stdout.strip().split("\n")

# Initialize CSV writer
csv_file = open('info.csv', 'w', newline='')
csv_writer = csv.writer(csv_file)
csv_writer.writerow(['Identifier', 'BioSample_Name'])

# Special case: handle files in base_directory
scp_base = [
    f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{args.base_directory}*.bam .",
    f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{args.base_directory}*.pbi .",
    f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{args.base_directory}*.xml ."
]
for cmd in scp_base:
    subprocess.run(cmd, shell=True)

# Iterate over directories
for dir_path in directories:
    # SCP commands
    scp_bam = f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{dir_path}*bam ."
    scp_bam_pbi = f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{dir_path}*bam.pbi ."
    scp_xml = f"scp -i {ssh_key_path} zachvanwinkle@136.165.158.10:{dir_path}*.xml ."
    
    # Execute SCP commands
    subprocess.run(scp_bam, shell=True)
    subprocess.run(scp_bam_pbi, shell=True)
    subprocess.run(scp_xml, shell=True)

    # Parse XML files
    for xml_file in dir_path.split('/')[-2:]:
        try:
            tree = ET.parse(xml_file)
            root = tree.getroot()
            bio_sample_name = root.find(".//pbsample:BioSample", namespaces={'pbsample': 'some_namespace'}).get('Name')
        except Exception as e:
            bio_sample_name = "N/A"
        
        # Write to CSV
        csv_writer.writerow([dir_path.split('/')[-2], bio_sample_name])

# Close CSV file
csv_file.close()
