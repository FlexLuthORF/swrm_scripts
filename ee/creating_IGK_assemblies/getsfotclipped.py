#!/bin/env python
import sys
import pysam

samfile = pysam.AlignmentFile(sys.argv[2])
bedfile = sys.argv[1]

def get_regions(filepath):
    regions = []
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            line = line.split("\t")
            line[1]=int(line[1])
            line[2]=int(line[2])
            regions.append(line)
    return regions
            
regions = get_regions(bedfile)

sumofSC = 0
counts = 0
bases = 0
frac = None
r_start = None
r_end = None

for region in regions:
    chrom = region[0]
    start = region[1]
    end = region[2]
    if r_start == None:
        r_start = start
    r_end = end
    for read in samfile.fetch(chrom, start, end):
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        counts = counts + 1
        sumofSC = read.get_cigar_stats()[0][4] + sumofSC
        bases = read.query_length + bases

if bases != 0:
    frac = sumofSC / bases

output = [chrom, r_start, r_end, bases, sumofSC, frac, counts]
output = map(str,output)
print ("\t".join(output))
        
# sc_bases = get_cigar_stats(bedfile, sum(samfile.get_cigar_stats()[4])

