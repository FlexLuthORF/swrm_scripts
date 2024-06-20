import os

# Function to copy files and rename them
def copy_and_rename(file_list_path):
    with open(file_list_path, 'r') as f:
        for line in f:
            sample_id, file_path = line.strip().split()
            file_name = os.path.basename(file_path)
            new_file_name = f"{sample_id}_{file_name}"
            os.system(f"cp {file_path} {new_file_name}")

# Path to the .txt file containing the file locations and sample IDs
file_list_path = "file_list.txt"

copy_and_rename(file_list_path)
