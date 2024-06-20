import os
import sqlite3
import pandas as pd
from pathlib import Path

# Function to create and connect to the SQLite database
def connect_db(db_name):
    conn = sqlite3.connect(db_name)
    return conn

# Function to setup the SQLite table
def setup_table(conn):
    cur = conn.cursor()
    cur.execute('''
        CREATE TABLE IF NOT EXISTS samples (
            unique_id TEXT PRIMARY KEY,
            sample_id TEXT,
            contigs_to_ref_phased_sorted_bam TEXT,
            ccs_to_ref_phased_sorted_bam TEXT,
            contigs_fasta TEXT,
            assembly_to_ref_bam TEXT,
            chr2_csv TEXT,
            chr2_imported_csv TEXT,
            chr22_csv TEXT,
            chr22_imported_csv TEXT,
            igh_csv TEXT,
            igh_imported_csv TEXT,
            ighc_csv TEXT,
            ighc_imported_csv TEXT,
            trb_csv TEXT,
            trb_imported_csv TEXT,
            ancestry_folder TEXT,
            inferred_ancestry TEXT,
            ancestry_likelihood REAL,
            aims_used TEXT
        )
    ''')
    conn.commit()

# Function to extract ancestry data
def extract_ancestry_data(file_path):
    if not Path(file_path).exists():
        print(f"File not found: {file_path}")
        return None, None, None
    try:
        with open(file_path, 'r') as file:
            content = file.readlines()
            inferred_ancestry = [line.split(': ')[1].strip() for line in content if "Inferred ancestry:" in line][0]
            ancestry_likelihood = float([line.split(': ')[1].strip() for line in content if inferred_ancestry + ':' in line][0])
            aims_used = [line.split(': ')[1].strip() for line in content if "Number of AIMs used:" in line][0]
            return inferred_ancestry, ancestry_likelihood, aims_used
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return None, None, None

# Function to search for files and extract data
def search_files_and_insert_data(sample_id, root_dir, conn):
    unique_id = f"{os.path.basename(root_dir)}_{sample_id}"
    file_paths = {
        'contigs_to_ref_phased_sorted_bam': f'{root_dir}/run_igenotyper/{sample_id}/alignments/contigs_to_ref_phased.sorted.bam',
        'ccs_to_ref_phased_sorted_bam': f'{root_dir}/run_igenotyper/{sample_id}/alignments/ccs_to_ref_phased.sorted.bam',
        'contigs_fasta': f'{root_dir}/run_igenotyper/{sample_id}/assembly/contigs.fasta',
        'assembly_to_ref_bam': f'{root_dir}/run_igenotyper/{sample_id}/assembly_to_ref.bam',
        # Initialize all potential keys with None to ensure correct bindings number
        'chr2_csv': None, 'chr2_imported_csv': None, 'chr22_csv': None, 'chr22_imported_csv': None,
        'igh_csv': None, 'igh_imported_csv': None, 'ighc_csv': None, 'ighc_imported_csv': None,
        'trb_csv': None, 'trb_imported_csv': None
    }

    annotations = ['chr2', 'chr22', 'igh', 'ighc', 'trb']
    for anno in annotations:
        csv_path = f'{root_dir}/run_igenotyper/{sample_id}/annotations/{sample_id}/{anno}/{sample_id}_make_gene_file.csv'
        imported_csv_path = f'{root_dir}/run_igenotyper/{sample_id}/annotations/{sample_id}/{anno}/{sample_id}_make_gene_file_imported.csv'
        if Path(csv_path).exists():
            file_paths[f'{anno}_csv'] = csv_path
        if Path(imported_csv_path).exists():
            file_paths[f'{anno}_imported_csv'] = imported_csv_path

    ancestry_folder = f'{root_dir}/run_igenotyper/{sample_id}/{sample_id}_Ancestry'
    ancestry_result_file = f'{ancestry_folder}/Sample_ancestry_result.txt'
    inferred_ancestry, ancestry_likelihood, aims_used = extract_ancestry_data(ancestry_result_file)

    # Insert data into SQLite
    cur = conn.cursor()
    cur.execute('''
        INSERT INTO samples (
            unique_id, sample_id, contigs_to_ref_phased_sorted_bam, ccs_to_ref_phased_sorted_bam, contigs_fasta,
            assembly_to_ref_bam, chr2_csv, chr2_imported_csv, chr22_csv, chr22_imported_csv, igh_csv, igh_imported_csv,
            ighc_csv, ighc_imported_csv, trb_csv, trb_imported_csv, ancestry_folder, inferred_ancestry,
            ancestry_likelihood, aims_used
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', [unique_id, sample_id] + [file_paths[key] for key in sorted(file_paths.keys())] + [ancestry_folder, inferred_ancestry, ancestry_likelihood, aims_used])
    conn.commit()


# Main function to run the script
def main(fofn_path, root_dir):
    conn = connect_db('sample_data.db')
    setup_table(conn)
    df = pd.read_csv(fofn_path, delimiter='\t', header=None, names=['sample_id'], usecols=[0])
    for sample_id in df['sample_id'].str.strip():
        search_files_and_insert_data(sample_id, root_dir, conn)
    conn.close()

# Specify the path to your file of names and the root directory
fofn_path = 'fofn.tsv'
root_dir = '/home/zmvanw01/projects/EF/240411'
main(fofn_path, root_dir)
