#!/bin/env python
import sys
#import vcf
import pandas as pd
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

vcffn = sys.argv[1]
samples = sys.argv[2]
matrixfn = sys.argv[3]
type_ = sys.argv[4]

def read_samples(fn):
    samples = []
    with open(fn,'r') as fh:
        for line in fh:
            samples.append(line.rstrip())
    return samples

def transform_genotype(genotype,type_):
    t = "NA"
    genotype = genotype.replace("|","/")
    if type_ == "biallelic":
        function = {
            "0/0": "0",
            "0/1": "1",
            "1/0": "1",
            "1/1": "2"
        }
    elif type_ == "hemi":
        function = {
            "./0": "0",
            "./1": "1",
            "1/.": "1",
            "0/.": "0"
        }
    elif type_ == "both":
        function = {
            "0/0": "0",
            "0/1": "1",
            "1/0": "1",
            "1/1": "2",
            "./0": "3",
            "0/.": "3",
            "./1": "4",
            "1/.": "4"
        }
    else:
        sys.exit("Wrong - 1")
    if genotype in function:
        t = function[genotype]
    # else:
    #     sys.exit(genotype)
    # if t == "NA":
    #     print(genotype)
    return t

def vcf_genotypes(fn, samples, type_):
    header = None
    sample_order = None
    sample_genotypes = [["chrom_pos"] + samples]
    with open(fn, 'r') as fh:
        for line in fh:
            if "##" in line:
                continue
            line = line.rstrip().split('\t')
            if "#" in line[0]:
                header = line
                sample_order = line[9:]
                continue
            gts = line[9:]
            sample_gts = {}
            for s, g in zip(sample_order, gts):
                sample_gts[s] = g.split(":")[0]
            t_gts = []
            for s in samples:
                t_gts.append(transform_genotype(sample_gts[s], type_))
            sample_genotype_counter = Counter(t_gts)
            # Modified part: concatenating chromosome and position
            chrom_pos = "{}_{}".format(line[0], line[1])
            pos_gts = [chrom_pos]
            output_gts = set()
            for s in samples:
                t_gt = transform_genotype(sample_gts[s],type_)
                if sample_genotype_counter[t_gt] < 2:
                    t_gt = "NA"
                pos_gts.append(t_gt)
                output_gts.add(t_gt)
            output_gts_count = len(output_gts)
            if "NA" in output_gts:
                output_gts_count = output_gts_count - 1
            if output_gts_count < 2: # 3
                continue
            sample_genotypes.append(pos_gts)
    return sample_genotypes

def process_sample(sample, vcffn, type_):
    genotypes = vcf_genotypes(vcffn, [sample], type_)
    return genotypes[1:]  # Exclude the header row

def main():
    if len(sys.argv) != 5:
        print("Usage: python script.py <vcf_file> <samples_file> <matrix_file> <type>")
        sys.exit(1)

    vcffn = sys.argv[1]
    samples_file = sys.argv[2]
    matrixfn = sys.argv[3]
    type_ = sys.argv[4]

    samples = read_samples(samples_file)

    with ThreadPoolExecutor(max_workers=12) as executor:
        futures = [executor.submit(process_sample, sample, vcffn, type_) for sample in samples]

        genotypes = [["chrom_pos"] + samples]
        for future in futures:
            sample_genotypes = future.result()
            for i in range(len(sample_genotypes)):
                if len(genotypes) <= i + 1:
                    genotypes.append(sample_genotypes[i])
                else:
                    genotypes[i + 1].extend(sample_genotypes[i][1:])

    with open(matrixfn, 'w') as fh:
        for row in genotypes:
            fh.write("%s\n" % "\t".join(row))

if __name__ == "__main__":
    main()