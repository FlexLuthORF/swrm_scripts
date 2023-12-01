import subprocess

# Define file paths
v_alleles_file = "path/to/V_alleles.fasta"  # Ungapped V alleles
d_alleles_file = "path/to/D_alleles.fasta"  # D alleles
j_alleles_file = "path/to/J_alleles.fasta"  # J alleles

# Gapping V alleles
# Define the gapped reference file and output file for gapped V alleles
gapped_ref_file = "path/to/gapped_ref.fasta"  # Reference for gapping
gapped_v_output = "path/to/gapped_V_alleles.fasta"  # Output for gapped V alleles
subprocess.run(["gap_sequences", v_alleles_file, gapped_ref_file, gapped_v_output])

# Create NDM file
# The ndm file is created from the gapped V alleles
ndm_output = "path/to/ndm_file.ndm"
subprocess.run(["make_igblast_ndm", gapped_v_output, "VH", ndm_output])  # Adjust "VH" if needed

# Create AUX file
# The aux file is created from J alleles
aux_output = "path/to/aux_file.aux"
subprocess.run(["annotate_j", j_alleles_file, aux_output])

# Create custom BLAST databases
subprocess.run(["makeblastdb", "-in", gapped_v_output, "-dbtype", "nucl"])
subprocess.run(["makeblastdb", "-in", d_alleles_file, "-dbtype", "nucl"])
subprocess.run(["makeblastdb", "-in", j_alleles_file, "-dbtype", "nucl"])

# IgBlast command with custom database
igblast_command = [
    "igblastn",
    "-germline_db_V", gapped_v_output,
    "-germline_db_D", d_alleles_file,
    "-germline_db_J", j_alleles_file,
    "--auxiliary_data", aux_output,
    "--custom_internal_data", ndm_output, 
]

igblast_output = "igblast_output.txt"
changeo_output = "changeo_output.txt"


subprocess.run(["MakeDb.py", "-i", igblast_output, "-s", sample_seq_file, "-o", changeo_output])


