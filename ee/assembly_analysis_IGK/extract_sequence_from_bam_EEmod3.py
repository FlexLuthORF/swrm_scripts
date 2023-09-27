#!/bin/env python
import sys
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

coordsfn = sys.argv[1]
bamfile = sys.argv[2]
coords = []

def get_read_coords(read, start, end):
    read_start = 0    
    read_end = None
    for query_pos, ref_pos in read.get_aligned_pairs():
        if query_pos is None:
            continue
        if ref_pos is None:
            continue
        if ref_pos < start:
            read_start = query_pos
        if ref_pos <= end:
            read_end = query_pos
        else:
            break
    assert read_end is not None
    return (read_start, read_end)

with open(coordsfn, 'r') as fh:
    for line in fh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        gene_name = line[3]  # Extract gene name from the 4th column
        coords.append([chrom, start, end, gene_name])

samfile = pysam.AlignmentFile(bamfile)
for chrom, start, end, gene_name in coords:
    reverse_complement = False
    
    if start < 89607941 and gene_name not in ['IGKV4-1', 'IGKV5-2']:
        reverse_complement = True
    
    for read in samfile.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.reference_end < start or read.reference_start > end:
            continue
        read_start, read_end = get_read_coords(read, start, end)
        read_sequence = read.query_sequence[read_start:read_end]
        
        if reverse_complement:
            read_sequence = str(Seq(read_sequence).reverse_complement())
        
        read_name = read.query_name
        haplotype = read.get_tag("RG").split(":")[-1]  # Extract haplotype from @RG tag
        fasta_header = f">{read_name}_{chrom}:{start}-{end}_{haplotype}_{gene_name}"  # Construct FASTA header
        print(fasta_header)
        print(read_sequence)

