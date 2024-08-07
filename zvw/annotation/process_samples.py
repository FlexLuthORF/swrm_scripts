import os
import argparse
import concurrent.futures
import subprocess

def process_sample(sample_id, bam_path):
    loci_options = [
        ('IGH', 'igh'),
        ('IGH', 'ighc'),
        ('IGK', 'chr2'),
        ('IGL', 'chr22'),
        ('TRB', 'trb'),
        ('TRG', 'chr7'),
        ('TRD', 'chr14'),
        ('TRA', 'chr14')
    ]
    # Get the directory above the bam_path
    output_base_dir = os.path.dirname(os.path.dirname(bam_path))
    # Run make_gene_file.py for each locus option
    for locus, loci in loci_options:
        output_dir = f"{output_base_dir}/annotations/{sample_id}/{locus}"
        os.makedirs(output_dir, exist_ok=True)
        make_gene_file_output = f"{output_dir}/{sample_id}_make_gene_file.csv"
        command = f"python /home/zmvanw01/git_repos/VDJbase-Genomics/python/make_gene_file.py --sample {sample_id} {locus} {loci} \"+-\" /home/zmvanw01/git_repos/immune_receptor_genomics/current /home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta {make_gene_file_output} -b {bam_path}"
        subprocess.run(command, shell=True, check=True)

    # Run import_from_assemblies.py for each locus option and merge the results
    merged_outfile = f"{output_base_dir}/annotations/{sample_id}/merged_annotations.csv"
    with open(merged_outfile, 'w') as outfile:
        for locus, loci in loci_options:
            gene_file = f"{output_base_dir}/annotations/{sample_id}/{locus}/{sample_id}_make_gene_file.csv"
            import_from_assembly_output = f"{output_base_dir}/annotations/{sample_id}/{locus}/{sample_id}_make_gene_file_imported.csv"
            command = f"python /home/zmvanw01/git_repos/VDJbase-Genomics/python/import_from_assemblies.py {locus} {loci} \"+-\" {gene_file} /home/zmvanw01/git_repos/immune_receptor_genomics/current/reference.fasta /home/zmvanw01/git_repos/immune_receptor_genomics/current /home/zmvanw01/git_repos/VDJbase-Genomics/ref {import_from_assembly_output}"
            subprocess.run(command, shell=True, check=True)
            with open(import_from_assembly_output, 'r') as infile:
                outfile.write(infile.read())

def main(samples):
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
        futures = [executor.submit(process_sample, sample_id, bam_path) for sample_id, bam_path in samples]
        concurrent.futures.wait(futures)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process samples using make_gene_file.py and import_from_assemblies.py')
    parser.add_argument('samples', nargs='+', help='sample ID and BAM path pairs')
    args = parser.parse_args()
    samples = [(args.samples[i], args.samples[i+1]) for i in range(0, len(args.samples), 2)]
    main(samples)