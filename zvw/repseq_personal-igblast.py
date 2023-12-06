import subprocess
import pandas as pd
import os

os.environ['IGDATA'] = '/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref'

#make if statements for loci with d_alleles

v_alleles_file = "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref/fasta/human_ig_V_alleles.fasta"  # Ungapped V alleles
#d_alleles_file = "path/to/D_alleles.fasta"  # D alleles
j_alleles_file = "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref/fasta/human_ig_J_alleles.fasta"  # J alleles
sample_seq_file = "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/S5-final_total_collapse-unique_atleast-2_reheader.fasta"  # Test sequences

# cp aux and optional from ~/share/igblast


# Create custom BLAST databases
#subprocess.run(["makeblastdb", "-in", v_alleles_file,
#                 "-dbtype", "nucl",
#                 "-parse_seqids",
#                 "-out", "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref/database/human_ig_V"])
#subprocess.run(["makeblastdb", "-in", d_alleles_file, "-dbtype", "nucl"])
#subprocess.run(["makeblastdb", "-in", j_alleles_file, 
#                "-dbtype", "nucl",
#                "-parse_seqids",
#                "-out", "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref/database/human_ig_J"])


igblast_output = "igblast_output.fmt7"
changeo_output = "changeo_output.tsv"

igblast_command = [
    "igblastn",
    "-germline_db_V", "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref/database/human_ig_V",
    "-germline_db_D", "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref/database/imgt_human_ig_d",
    "-germline_db_J", "/home/zmvanw01/projects/CW50/231122/presto/CW50_NS_pool1_220781702C/alleles/personal-ref/database/human_ig_J",
    "-auxiliary_data", aux_output,
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

