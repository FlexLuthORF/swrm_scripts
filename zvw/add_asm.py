import sys
import os
import pandas as pd

def add_assembly_path(tsv_path, project_folder, output_path):
    # Load the TSV file
    df = pd.read_csv(tsv_path, sep='\t', header=None)
    
    # Define the assembly path using the sampleId and check its existence
    def get_assembly_path(sampleId):
        path = f"/home/zmvanw01/projects/{project_folder}/run_hifiasm/{sampleId}/merged_bam/alg_asm20_to_ref_with_secondarySeq/{sampleId}.sorted.bam"
        return path if os.path.exists(path) else "Path does not exist"

    df[2] = df[0].apply(get_assembly_path)
    
    # Write the modified DataFrame to a new TSV file
    df.to_csv(output_path, sep='\t', header=False, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_tsv_path> <project_folder> <output_tsv_path>")
        sys.exit(1)
    
    tsv_path = sys.argv[1]
    project_folder = sys.argv[2]
    output_path = sys.argv[3]
    
    add_assembly_path(tsv_path, project_folder, output_path)
