import subprocess
import pandas as pd
import os
import argparse
import glob

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('sampleID', type=str, help='Sample ID')
args = parser.parse_args()

sampleID = args.sampleID

for chain in ('IGK','IGL'):
    
    database_dir = f'/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database'
    fasta_dir = f'/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/fasta'
    igdata_path=f'/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}'
    os.environ['IGDATA'] = igdata_path
    os.makedirs(igdata_path, exist_ok=True)
    os.makedirs(database_dir, exist_ok=True)
    os.makedirs(fasta_dir, exist_ok=True)
    # cp aux and optional from ~/share/igblast
    subprocess.run(["cp", "-r", "/home/zmvanw01/share/igblast/internal_data", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}"])
    subprocess.run(["cp", "-r", "/home/zmvanw01/share/igblast/optional_file", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}"])
    d_allele_files = glob.glob('/home/zmvanw01/share/igblast/database/imgt_human_ig_d.*')
    for file in d_allele_files:
        subprocess.run(["cp", file, f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database"])



    subprocess.run(["python",
                    "/home/zmvanw01/git_repos/swrm_scripts/zvw/fasta-from-annotations.py", 
                   f"/home/egenge01/projects/CW50/{chain}_alleles/{sampleID}/annotations.csv", 
                   f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}"])
    
    v_alleles_file = f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/fasta/human_ig_V_alleles.fasta"  # Ungapped V alleles
    d_alleles_file = f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/fasta/human_ig_D_alleles.fasta"  # D alleles
    j_alleles_file = f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/fasta/human_ig_J_alleles.fasta"  # J alleles
    sample_seq_file = f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/S5-final_total_collapse-unique_atleast-2_reheader.fasta"  # Test sequences

    # Create custom BLAST databases
    subprocess.run(["makeblastdb", "-in", v_alleles_file,
                    "-dbtype", "nucl",
                    "-parse_seqids",
                    "-out", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database/human_ig_V"])
   # if chain != "IGH":
    #    subprocess.run(["makeblastdb", "-in", d_alleles_file, "-dbtype", "nucl" "-parse_seqids", "-out", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database/human_ig_D"])
    subprocess.run(["makeblastdb", "-in", j_alleles_file, 
                    "-dbtype", "nucl",
                    "-parse_seqids",
                    "-out", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database/human_ig_J"])
    
    changeo_folder=f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/changeo"
    os.makedirs(igdata_path, exist_ok=True)
    igblast_output = os.path.join(changeo_folder, f"igblast_output_{chain}.fmt7")

    igblast_command = [
        "igblastn",
        "-germline_db_V", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database/human_ig_V",
        "-germline_db_D", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database/imgt_human_ig_d",
        "-germline_db_J", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/database/human_ig_J",
        "-auxiliary_data", f"/home/zmvanw01/projects/CW50/231122/presto/{sampleID}/alleles/personal-ref/{chain}/optional_file/human_gl.aux",
        #"-custom_internal_data", ndm_output, 
        "-query", sample_seq_file,
        "-outfmt", "19",
        "-out", igblast_output,
        #"-num_threads", "4",
        "-organism", "human",
        "-ig_seqtype", "Ig",
        "-domain_system", "imgt"
    ]

    subprocess.run(igblast_command)

#now handle IGH
igdata_path='/home/zmvanw01/share/IGH'
os.environ['IGDATA'] = igdata_path
igblast_output = os.path.join(changeo_folder, f"igblast_output_H.fmt7")
igblast_command = [
        "igblastn",
        "-germline_db_V", "/home/zmvanw01/share/IGH/database/human_ig_V",
        "-germline_db_D", "/home/zmvanw01/share/IGH/database/imgt_human_ig_d",
        "-germline_db_J", "/home/zmvanw01/share/IGH/database/human_ig_J",
        "-auxiliary_data", "/home/zmvanw01/share/IGH/optional_file/human_gl.aux",
        #"-custom_internal_data", ndm_output, 
        "-query", sample_seq_file,
        "-outfmt", "19",
        "-out", igblast_output,
        #"-num_threads", "4",
        "-organism", "human",
        "-ig_seqtype", "Ig",
        "-domain_system", "imgt"
    ]
