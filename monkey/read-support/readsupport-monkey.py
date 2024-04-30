import os
import subprocess
import sys

def run_map_ccs_to_pers(scratch, sample, ccs_bam, assemblies_fasta, monkey_mask_ref):
    outd = os.path.join(scratch, "read_support", sample)
    os.makedirs(os.path.join(outd, "ccs_to_pers"), exist_ok=True)

    # Convert PacBio HiFi reads to FASTA
    with open(os.path.join(outd, "ccs_to_pers", "reads.fasta"), "w") as reads_fasta:
        samtools_view = subprocess.Popen(["samtools", "view", ccs_bam], stdout=subprocess.PIPE)
        awk_process = subprocess.Popen(["awk", "{ print \">\"$1\"\\n\"$10 }"], stdin=samtools_view.stdout, stdout=reads_fasta)
        samtools_view.stdout.close()
        awk_process.communicate()

    subprocess.run(["samtools", "faidx", os.path.join(outd, "ccs_to_pers", "reads.fasta")], check=True)

    # Create personalized reference
    with open(os.path.join(outd, "ccs_to_pers", "contigs.fasta"), "w") as contigs_fasta:
        samtools_view = subprocess.Popen(["samtools", "view", "-F", "0x100", "-F", "0x800", assemblies_fasta], stdout=subprocess.PIPE)
        awk_process = subprocess.Popen(["awk", "{print \">\"$1\"\\n\"$10}"], stdin=samtools_view.stdout, stdout=contigs_fasta)
        samtools_view.stdout.close()
        awk_process.communicate()

    subprocess.run(["samtools", "faidx", os.path.join(outd, "ccs_to_pers", "contigs.fasta")], check=True)

    pers_contigs = os.path.join(outd, "ccs_to_pers", "contigs.fasta")

    with open(os.path.join(outd, "ccs_to_pers", "pers_ref.fasta"), "w") as pers_ref_fasta:
        subprocess.run(["cat", monkey_mask_ref, pers_contigs], stdout=pers_ref_fasta, check=True)

    subprocess.run(["samtools", "faidx", os.path.join(outd, "ccs_to_pers", "pers_ref.fasta")], check=True)

    with open(os.path.join(outd, "ccs_to_pers", "output.sam"), "w") as output_sam:
        subprocess.run(["minimap2", "-ax", "map-hifi", "--secondary=no", "-t", "10", "-L", os.path.join(outd, "ccs_to_pers", "pers_ref.fasta"), os.path.join(outd, "ccs_to_pers", "reads.fasta")], stdout=output_sam, check=True)

    with open(os.path.join(outd, "ccs_to_pers", "output.bam"), "w") as output_bam:
        subprocess.run(["samtools", "view", "-Sbh", os.path.join(outd, "ccs_to_pers", "output.sam")], stdout=output_bam, check=True)

    subprocess.run(["samtools", "sort", "-@", "10", os.path.join(outd, "ccs_to_pers", "output.bam"), "-o", os.path.join(outd, "ccs_to_pers", "output.sorted.bam")], check=True)
    subprocess.run(["samtools", "index", os.path.join(outd, "ccs_to_pers", "output.sorted.bam")], check=True)
    os.remove(os.path.join(outd, "ccs_to_pers", "output.sam"))

def get_read_support_vdj3(scratch, sample, igh_digger, igk_digger, igl_digger):
    base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")
    bam_file = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "output.sorted.bam")
    ref = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "pers_ref.fasta")

    # Create the base output directory if it doesn't exist
    os.makedirs(base_outd, exist_ok=True)

    if not os.path.exists(f"{bam_file}.bai"):
        subprocess.run(["samtools", "index", bam_file], check=True)

    for gene_type, import_out in [("igh", igh_digger), ("igk", igk_digger), ("igl", igl_digger)]:
        if os.path.exists(import_out):
            os.makedirs(os.path.join(base_outd, gene_type), exist_ok=True)

            tmp_file = os.path.join(scratch, "read_support", sample, "tmp", f"{os.path.basename(import_out)}_read_support.tmp")
            os.makedirs(os.path.dirname(tmp_file), exist_ok=True)

            with open(tmp_file, "w") as tmp_file_write:
                tmp_file_write.write("Total_Positions,Average_Coverage,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Percent_Accuracy,Positions_With_At_Least_10x_Coverage,Fully_Spanning_Reads,Fully_Spanning_Reads_100%_Match\n")

            with open(import_out, "r") as import_out_read:
                header = import_out_read.readline().strip()
                header_cols = header.split(",")
                contig_col = header_cols.index("contig")
                start_col = header_cols.index("start")
                end_col = header_cols.index("end")
                gene_col = header_cols.index("ASC_match")

            with open(import_out, "r") as import_out_read:
                next(import_out_read)  # Skip header
                for line in import_out_read:
                    line = line.strip().split(",")
                    contig = line[contig_col]
                    start = int(float(line[start_col]))
                    end = int(float(line[end_col]))
                    gene = line[gene_col]
                    region = f"{contig}:{start}-{end}"

                    contig_filename = contig.replace("/", "_")
                    tmp_bam = os.path.join(base_outd, gene_type, f"{contig_filename}_{start}_{end}.bam")
                    os.makedirs(os.path.dirname(tmp_bam), exist_ok=True)
                    subprocess.run(["samtools", "view", "-F", "0x100", "-F", "0x800", "-b", bam_file, "-o", tmp_bam, "-U", "/dev/null", region], check=True)
                    subprocess.run(["samtools", "index", tmp_bam], check=True)

                    with open(f"{tmp_file}_awk_out", "w") as awk_out:
                        samtools_mpileup = subprocess.Popen(["samtools", "mpileup", "-f", ref, "-r", region, tmp_bam], stdout=subprocess.PIPE)
                        awk_process = subprocess.Popen(["awk", "--csv", "-v", f"total_positions={end - start + 1}", "-v", f"sample={sample}", 'BEGIN { total_reads=0; mismatched_positions=0; matched_positions=0; positions_with_10x=0; mismatch_list=""; match_list=""; } { total_reads += length($5); mismatches = length(gensub(/[.,]/, "", "g", $5)); matches = length(gensub(/[^.,]/, "", "g", $5)); mismatch_list = (mismatch_list == "" ? mismatches : mismatch_list ":" mismatches); match_list = (match_list == "" ? matches : match_list ":" matches); coverage = length($5); if (coverage >= 10) { positions_with_10x++; } mismatch_rate = mismatches / coverage; match_rate = matches / coverage; if (mismatch_rate > 0.2) { mismatched_positions++; } if (match_rate > 0.8) { matched_positions++; } } END { avg_reads_per_position = (total_positions > 0) ? total_reads / total_positions : 0; percent_accuracy = (matched_positions / total_positions) * 100; print total_positions, avg_reads_per_position, mismatched_positions, matched_positions, mismatch_list, match_list, percent_accuracy, positions_with_10x; }'], stdin=samtools_mpileup.stdout, stdout=awk_out)
                        samtools_mpileup.stdout.close()
                        awk_process.communicate()

                    subprocess.run(["python", "/home/zmvanw01/git_repos/swrm_scripts/monkey/read-support/match_subsequences.py", tmp_bam, contig, str(start), str(end), gene, import_out], stdout=open(f"{tmp_file}_py_out", "w"), check=True)

                    print(f"gene is {gene}. awk_out is {open(f'{tmp_file}_awk_out').read().strip()}. py_out is {open(f'{tmp_file}_py_out').read().strip()}")

                    with open(tmp_file, "a") as tmp_file_append, open(f"{tmp_file}_awk_out", "r") as awk_out_read, open(f"{tmp_file}_py_out", "r") as py_out_read:
                        tmp_file_append.write(f"{awk_out_read.read().strip()},{py_out_read.read().strip()}\n")

                    os.remove(f"{tmp_file}_awk_out")
                    os.remove(f"{tmp_file}_py_out")
                    os.remove(tmp_bam)
                    os.remove(f"{tmp_bam}.bai")

            if os.path.exists(import_out) and os.path.exists(tmp_file):
                combined_file = os.path.join(scratch, "read_support", sample, "output", gene_type, f"{os.path.splitext(os.path.basename(import_out))[0]}_combined.csv")
                os.makedirs(os.path.dirname(combined_file), exist_ok=True)
                with open(combined_file, "w") as combined_file_write, open(import_out, "r") as import_out_read, open(tmp_file, "r") as tmp_file_read:
                    for line1, line2 in zip(import_out_read, tmp_file_read):
                        combined_file_write.write(f"{line1.strip()},{line2.strip()}\n")

if __name__ == "__main__":
    sample = sys.argv[1]
    assemblies_fasta = sys.argv[2]
    igh_digger = sys.argv[3]
    igk_digger = sys.argv[4]
    igl_digger = sys.argv[5]
    ccs_bam = sys.argv[6]
    monkey_mask_ref = "/home/zmvanw01/ref/monkey/Child_17thApril2024_bothhap_IG_masked_unique.fasta"
    scratch = os.getcwd()

    run_map_ccs_to_pers(scratch, sample, ccs_bam, assemblies_fasta, monkey_mask_ref)
    get_read_support_vdj3(scratch, sample, igh_digger, igk_digger, igl_digger)