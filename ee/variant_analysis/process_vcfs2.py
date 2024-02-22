#!/usr/bin/env python2
import sys
import os
import shutil
import subprocess

def read_data(filename):
    with open(filename, 'r') as file:
        headers = file.readline().strip().split('\t')
        data = [line.strip().split('\t') for line in file]
    return headers, data

def read_bed_file(bed_filename):
    bed_entries = {}
    with open(bed_filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            bed_entries[parts[3]] = parts[:3]  # Store chromosome, start, and end
    return bed_entries

def create_sample_bed_file(hemi_regions, bed_entries, sample_id, outd):
    sample_bed_filename = os.path.join(outd, '{}_hemi.bed'.format(sample_id))
    with open(sample_bed_filename, 'w') as bed_file:
        for region in hemi_regions:
            bed_file.write('\t'.join(bed_entries[region]) + '\t' + region + '\n')
    return sample_bed_filename

def process_samples(headers, data, bed_entries, outd):
    for row in data:
        sample_id = row[0]
        hemi_regions = [headers[i] for i, genotype in enumerate(row[1:], 1) if genotype == '0/1' and headers[i] in bed_entries]

        output_vcf_dir = os.path.join(outd, sample_id)
        if not os.path.exists(output_vcf_dir):
            os.makedirs(output_vcf_dir)
        output_vcf = os.path.join(output_vcf_dir, '{}.vcf'.format(sample_id))
        input_vcf = '/home/egenge01/projects/CW50/vcfs/{}/{}_call_norm.vcf'.format(sample_id, sample_id)

        if hemi_regions:
            sample_bed_filename = create_sample_bed_file(hemi_regions, bed_entries, sample_id, outd)
            cmd = ["python", "/home/egenge01/bioinformatics/python/change_genotypes.py", 
                   input_vcf, sample_bed_filename, output_vcf]
            subprocess.call(cmd)
        else:
            # Copy the original VCF file to the output directory as is
            shutil.copy(input_vcf, output_vcf)

def main(genotypes_filename, bed_filename):
    outd = os.path.join(os.getcwd(), 'vcf_change_to_hemi')
    headers, data = read_data(genotypes_filename)
    bed_entries = read_bed_file(bed_filename)
    process_samples(headers, data, bed_entries, outd)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage: python process_samples.py <genotypes_filename> <bed_filename>"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
