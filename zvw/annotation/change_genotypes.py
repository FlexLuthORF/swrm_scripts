#!/usr/bin/env python3
import vcf
import sys

from collections import namedtuple

def create_new_record(record, samples_updated, samples_indexes_updated):
    new_record = vcf.model._Record(record.CHROM,
                                  record.POS,
                                  record.ID,
                                  record.REF,
                                  record.ALT,
                                  record.QUAL,
                                  record.FILTER,
                                  record.INFO,
                                  record.FORMAT,
                                  samples_updated,
                                  samples_indexes_updated)
    return new_record

def read_bedfile(bedfile):
    coords = []
    with open(bedfile, 'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            coords.append([chrom, start, end])
    return coords

def main():
    if len(sys.argv) != 4:
        print("Usage: {} <input_vcf> <heterozygous_deletions> <output_path>".format(sys.argv[0]))
        sys.exit(1)

    input_vcf = sys.argv[1]
    bedfile = sys.argv[2]
    output_path = sys.argv[3]

    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    vcf_writer = vcf.Writer(open(output_path, 'w'), vcf_reader)

    coords = read_bedfile(bedfile)

    CallData = namedtuple('CallData', ['GT'])
    samples = vcf_reader.samples

    for record in vcf_reader:
        samples_updated = []
        samples_indexes_updated = {}
        new_record = create_new_record(record, samples_updated, samples_indexes_updated)

        for index, sample in enumerate(samples):
            change_snp = False
            for chrom, start, end in coords:
                if int(record.POS) < start:
                    continue
                if int(record.POS) > end:
                    continue
                change_snp = True
                break

            if change_snp:
                current_genotype = record.genotype(sample)['GT']
                if "/" in current_genotype:
                    alleles = current_genotype.split("/")
                else:
                    alleles = current_genotype.split("|")
                if "." in alleles:
                    samples_updated.append(record.samples[index])
                    continue
                alleles = map(int, alleles)
                new_genotype = "./%s" % max(alleles)
                c = CallData(GT=new_genotype)
                new_call = vcf.model._Call(new_record, sample, c)
                samples_updated.append(new_call)
            else:
                samples_updated.append(record.samples[index])
            samples_indexes_updated[sample] = index
        new_record._sample_indexes = samples_indexes_updated
        new_record.samples = samples_updated
        vcf_writer.write_record(new_record)

if __name__ == "__main__":
    main()
