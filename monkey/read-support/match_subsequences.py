#!/usr/bin/env python3

import csv
import sys
import pysam

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in seq[::-1]])

# Function to extract the correct sequence to compare against
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

# Function to count matching reads
def count_matching_reads(bamfile, chrom, start, end, sequence):
    full_span_count = 0
    perfect_match_count = 0
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for read in samfile.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.reference_start <= start and read.reference_end >= end:
            full_span_count += 1
            read_seq = read.query_sequence
            if sequence in read_seq or reverse_complement(sequence) in read_seq:
                perfect_match_count += 1
    return full_span_count, perfect_match_count

# Main script execution
if __name__ == '__main__':
    bamfile = sys.argv[1]
    contig = sys.argv[2]
    start = int(sys.argv[3])
    end = int(sys.argv[4])
    gene_key = sys.argv[5]
    import_csv = sys.argv[6]  # The CSV file with the sequences

    sequence = get_sequence_from_csv(import_csv, gene_key, contig, start, end)
    if sequence:
        full_span_count, perfect_match_count = count_matching_reads(bamfile, contig, start, end, sequence)
        print(f"{full_span_count},{perfect_match_count}")
    else:
        print(f"No sequence found for gene {gene_key} in region {contig}:{start}-{end}", file=sys.stderr)