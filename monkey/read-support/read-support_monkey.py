#!/usr/bin/env python3

import os
import sys
import csv
import subprocess
import tempfile

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])

def get_sequence_from_csv(import_csv, gene_key, contig, start, end):
    sequence = ""
    reverse_comp = False
    with open(import_csv, mode='r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            csv_start = int(float(row['start']))
            csv_end = int(float(row['end']))
            #print(f"gene is:{gene_key}. contig is: {contig}. start is: {start}. end is {end}.")
            if row['ASC_match'] == gene_key and row['contig'] == contig and csv_start == start and csv_end == end:
                if 'sense' in row and row['sense'] == '-':
                    reverse_comp = True

                sequence = row['seq']
                if reverse_comp:
                    sequence = reverse_complement(sequence)
                break

    return sequence

def count_matching_reads(bamfile, chrom, start, end, sequence):
    cmd = f"samtools view -F 0x100 -F 0x800 {bamfile} {chrom}:{start}-{end}"
    output = subprocess.check_output(cmd, shell=True, universal_newlines=True)

    full_span_count = 0
    perfect_match_count = 0
    for line in output.split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        read_start = int(fields[3])
        read_end = read_start + len(fields[9])
        if read_start <= start and read_end >= end:
            full_span_count += 1
            read_seq = fields[9]
            if sequence in read_seq or reverse_complement(sequence) in read_seq:
                perfect_match_count += 1

    return full_span_count, perfect_match_count

def run_make_ref_masked(monkey_mask_ref):
    subprocess.run(f"samtools faidx {monkey_mask_ref}", shell=True)

def run_map_ccs_to_pers(sample, scratch, ccs_bam, assemblies_fasta, monkey_mask_ref):
    outd = os.path.join(scratch, "read_support", sample)
    os.makedirs(os.path.join(outd, "ccs_to_pers"), exist_ok=True)

    with open(os.path.join(outd, "ccs_to_pers", "reads.fasta"), "w") as outfile:
        cmd = f"samtools view {ccs_bam} | awk '{{ print \">\"$1\"\\n\"$10 }}'"
        subprocess.run(cmd, shell=True, stdout=outfile)

    subprocess.run(f"samtools faidx {os.path.join(outd, 'ccs_to_pers', 'reads.fasta')}", shell=True)

    with open(os.path.join(outd, "ccs_to_pers", "contigs.fasta"), "w") as outfile:
        cmd = f"samtools view -F 0x100 -F 0x800 {assemblies_fasta} | awk '{{print \">\"$1\"\\n\"$10}}'"
        subprocess.run(cmd, shell=True, stdout=outfile)

    subprocess.run(f"samtools faidx {os.path.join(outd, 'ccs_to_pers', 'contigs.fasta')}", shell=True)

    pers_contigs = os.path.join(outd, "ccs_to_pers", "contigs.fasta")

    with open(os.path.join(outd, "ccs_to_pers", "pers_ref.fasta"), "w") as outfile:
        with open(monkey_mask_ref, "r") as infile:
            outfile.write(infile.read())
        with open(pers_contigs, "r") as infile:
            outfile.write(infile.read())

    # Uniqify the headers in pers_ref.fasta using awk
    #awk '/^>/ {++counter; print ">seq_" counter; next} {print}' ${outd}/ccs_to_pers/pers_ref.fasta > ${outd}/ccs_to_pers/pers_ref_uniqified.fasta
    #mv ${outd}/ccs_to_pers/pers_ref_uniqified.fasta ${outd}/ccs_to_pers/pers_ref.fasta

    #echo "Uniqified FASTA file: ${outd}/ccs_to_pers/pers_ref.fasta"

    subprocess.run(f"samtools faidx {os.path.join(outd, 'ccs_to_pers', 'pers_ref.fasta')}", shell=True)

    cmd = f"minimap2 -ax map-hifi --secondary=no -t 10 -L {os.path.join(outd, 'ccs_to_pers', 'pers_ref.fasta')} {os.path.join(outd, 'ccs_to_pers', 'reads.fasta')} | samtools view -Sbh - | samtools sort -@ 10 -o {os.path.join(outd, 'ccs_to_pers', 'output.sorted.bam')} && samtools index {os.path.join(outd, 'ccs_to_pers', 'output.sorted.bam')}"
    subprocess.run(cmd, shell=True)

    #rm -f ${outd}/ccs_to_pers/output.sam

def get_read_support_vdj3(sample, scratch, igh_digger, igk_digger, igl_digger):
    base_outd = os.path.join(scratch, "read_support", sample, "imported_genes")
    bam_file = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "output.sorted.bam")
    ref = os.path.join(scratch, "read_support", sample, "ccs_to_pers", "pers_ref.fasta")

    os.makedirs(base_outd, exist_ok=True)

    if not os.path.exists(f"{bam_file}.bai"):
        subprocess.run(f"samtools index {bam_file}", shell=True)

    gene_types = {
        "igh": igh_digger,
        "igk": igk_digger,
        "igl": igl_digger
    }

    for gene_type, import_out in gene_types.items():
        os.makedirs(os.path.join(base_outd, gene_type), exist_ok=True)

        if os.path.exists(import_out):
            with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
                tmp_file.write("Total_Positions,Average_Coverage,Mismatched_Positions,Matched_Positions,Position_Mismatches,Position_Matches,Percent_Accuracy,Positions_With_At_Least_10x_Coverage,Fully_Spanning_Reads,Fully_Spanning_Reads_100%_Match\n")
                tmp_file_path = tmp_file.name
            print(f"Created temporary file: {tmp_file_path}")

            header_cols = None
            with open(import_out, "r") as infile:
                header = next(infile).strip()
                header_cols = header.split(",")
                contig_col = header_cols.index("contig")
                start_col = header_cols.index("start")
                end_col = header_cols.index("end")
                gene_col = header_cols.index("ASC_match")

            with open(import_out, "r") as infile, tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_counts:
                next(infile)  # Skip header
                for line in infile:
                    fields = line.strip().split(",")
                    contig = fields[contig_col]
                    start = int(float(fields[start_col]))
                    end = int(float(fields[end_col]))
                    gene = fields[gene_col]

                    if start > end:
                        print(f"Skipping invalid region: {contig}:{start}-{end}", file=sys.stderr)
                        continue

                    region = f"{contig}:{start}-{end}"

                    contig_filename = contig.replace("/", "_")
                    tmp_bam = os.path.join(base_outd, gene_type, f"{contig_filename}_{start}_{end}.bam")
                    os.makedirs(os.path.dirname(tmp_bam), exist_ok=True)

                    cmd = f"samtools view -F 0x100 -F 0x800 -b {bam_file} -o {tmp_bam} -U /dev/null {region}"
                    subprocess.run(cmd, shell=True)
                    subprocess.run(f"samtools index {tmp_bam}", shell=True)

                    cmd = f"samtools mpileup -f {ref} -r {region} {tmp_bam}"
                    try:
                        output = subprocess.check_output(cmd, shell=True, universal_newlines=True, stderr=subprocess.DEVNULL)
                    except subprocess.CalledProcessError as e:
                        print(f"Error running samtools mpileup for region {region}: {e}", file=sys.stderr)
                        continue

                    total_positions = end - start + 1
                    total_reads = 0
                    mismatched_positions = 0
                    matched_positions = 0
                    positions_with_10x = 0
                    mismatch_list = []
                    match_list = []

                    for pileup_line in output.split("\n"):
                        if not pileup_line:
                            continue
                        fields = pileup_line.split("\t")
                        coverage = len(fields[4])
                        total_reads += coverage
                        mismatches = len(fields[4].replace(".", "").replace(",", ""))
                        matches = len(fields[4]) - mismatches

                        mismatch_list.append(str(mismatches))
                        match_list.append(str(matches))

                        if coverage >= 10:
                            positions_with_10x += 1

                        mismatch_rate = mismatches / coverage
                        match_rate = matches / coverage

                        if mismatch_rate > 0.2:
                            mismatched_positions += 1

                        if match_rate > 0.8:
                            matched_positions += 1

                    mismatch_list_str = ":".join(mismatch_list)
                    match_list_str = ":".join(match_list)

                    avg_reads_per_position = total_reads / total_positions if total_positions > 0 else 0
                    percent_accuracy = (matched_positions / total_positions) * 100

                    sequence = get_sequence_from_csv(import_out, gene, contig, start, end)
                    fully_spanning_reads, fully_spanning_reads_100_match = count_matching_reads(tmp_bam, contig, start, end, sequence)

                    tmp_counts.write(f"{total_positions},{avg_reads_per_position},{mismatched_positions},{matched_positions},{mismatch_list_str},{match_list_str},{percent_accuracy},{positions_with_10x},{fully_spanning_reads},{fully_spanning_reads_100_match}\n")

                    # python /home/zmvanw01/git_repos/swrm_scripts/monkey/read-support/match_subsequences.py "$tmp_bam" "$contig" "$start" "$end" "$gene" "$import_out" > "${tmp_file}_py_out"

                tmp_counts_path = tmp_counts.name
            print(f"Created temporary counts file: {tmp_counts_path}")

            if os.path.exists(import_out) and os.path.exists(tmp_file_path):
                combined_file = os.path.join(scratch, "read_support", sample, "output", gene_type, f"{os.path.basename(import_out).rstrip('.csv')}_combined.csv")
                final_output = os.path.join(scratch, "read_support", sample, "output", gene_type, f"{os.path.basename(import_out).rstrip('.csv')}_with_read_support.csv")
                os.makedirs(os.path.dirname(combined_file), exist_ok=True)
                print(f"Creating combined file: {combined_file}")
                print(f"Creating final output file: {final_output}")

                with open(combined_file, "w") as outfile:
                    with open(import_out, "r") as infile1, open(tmp_file_path, "r") as infile2:
                        for line1, line2 in zip(infile1, infile2):
                            outfile.write(f"{line1.strip()},{line2}")

                #rm -f "$tmp_file"

                with open(combined_file, "r") as infile, open(final_output, "w") as outfile:
                    reader = csv.reader(infile)
                    writer = csv.writer(outfile)
                    header = next(reader)
                    writer.writerow(header)
                    for row in reader:
                        row[1] = sample
                        row[2] = sample
                        writer.writerow(row)

                #rm -f "$combined_file"

            os.remove(tmp_file_path)
            os.remove(tmp_counts_path)
            #rm "${tmp_file}_awk_out" "${tmp_file}_py_out"
            #rm "$tmp_bam" "${tmp_bam}.bai"

def main():
    sample = sys.argv[1]
    assemblies_fasta = sys.argv[2]
    igh_digger = sys.argv[3]
    igk_digger = sys.argv[4]
    igl_digger = sys.argv[5]
    ccs_bam = sys.argv[6]
    monkey_mask_ref = "/home/zmvanw01/ref/monkey/Child_17thApril2024_bothhap_IG_masked_unique.fasta"

    scratch = os.getcwd()

    #run_make_ref_masked(monkey_mask_ref)
    run_map_ccs_to_pers(sample, scratch, ccs_bam, assemblies_fasta, monkey_mask_ref)
    get_read_support_vdj3(sample, scratch, igh_digger, igk_digger, igl_digger)

if __name__ == "__main__":
    main()