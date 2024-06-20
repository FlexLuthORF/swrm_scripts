import os
import sqlite3
import pandas as pd

# Function to create and connect to the SQLite database
def connect_db(db_name):
    conn = sqlite3.connect(db_name)
    return conn

# Function to ensure new columns exist in the table
def ensure_columns(conn):
    cur = conn.cursor()
    # Adding new columns for hifiasm data, if they don't already exist
    columns_to_add = [
        'contigs_to_ref_hifiasm',
        'chr2_hifiasm', 'chr2_imported_hifiasm', 'chr22_hifiasm', 'chr22_imported_hifiasm',
        'igh_hifiasm', 'igh_imported_hifiasm', 'ighc_hifiasm', 'ighc_imported_hifiasm',
        'trb_hifiasm', 'trb_imported_hifiasm'
    ]
    for col in columns_to_add:
        try:
            cur.execute(f'ALTER TABLE samples ADD COLUMN {col} TEXT')
        except sqlite3.OperationalError:
            # Column already exists; ignore the error
            pass
    conn.commit()

# Function to search for hifiasm files and update data
def search_files_and_update_data(sample_id, root_dir, conn):
    unique_id = f"{os.path.basename(root_dir)}_{sample_id}"
    file_paths = {
        'contigs_to_ref_hifiasm': f'{root_dir}/run_hifiasm/{sample_id}/merged_bam/alg_asm20_to_ref_with_secondarySeq/{sample_id}.sorted.bam'
    }

    annotations = ['chr2', 'chr22', 'igh', 'ighc', 'trb']
    for anno in annotations:
        file_paths[f'{anno}_hifiasm'] = f'{root_dir}/run_hifiasm/{sample_id}/merged_bam/annotations/{sample_id}/{anno}/{sample_id}_make_gene_file.csv'
        file_paths[f'{anno}_imported_hifiasm'] = f'{root_dir}/run_hifiasm/{sample_id}/merged_bam/annotations/{sample_id}/{anno}/{sample_id}_make_gene_file_imported.csv'

    # Update data into SQLite
    cur = conn.cursor()
    update_columns = ', '.join([f'{k} = ?' for k in file_paths.keys()])
    update_values = list(file_paths.values())
    update_values.append(unique_id)
    cur.execute(f'''
        UPDATE samples SET {update_columns} WHERE unique_id = ?
    ''', update_values)
    conn.commit()

# Main function to run the script
def main(fofn_path, root_dir, db_name):
    conn = connect_db(db_name)
    ensure_columns(conn)
    df = pd.read_csv(fofn_path, delimiter='\t', header=None, names=['sample_id'], usecols=[0])
    for sample_id in df['sample_id'].str.strip():
        search_files_and_update_data(sample_id, root_dir, conn)
    conn.close()

# Specify the path to your file of names and the root directory
fofn_path = 'fofn.tsv'
root_dir = '/home/zmvanw01/projects/EF/240411'
db_name = '240411.db'
main(fofn_path, root_dir, db_name)