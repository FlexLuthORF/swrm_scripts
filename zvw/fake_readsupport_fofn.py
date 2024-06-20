import sys

def process_fofn(input_file, output_file):
    fake_paths = {
        "chr2_gene": "/home/zmvanw01/projects/EF/240505/run_hifiasm/101-1027/merged_bam/annotations/101-1027/chr2/101-1027_make_gene_file.csv",
        "chr2_import": "/home/zmvanw01/projects/EF/240505/run_hifiasm/101-1027/merged_bam/annotations/101-1027/chr2/101-1027_make_gene_file_imported.csv",
        "chr22_gene": "/home/zmvanw01/projects/EF/240505/run_hifiasm/101-1027/merged_bam/annotations/101-1027/chr22/101-1027_make_gene_file.csv",
        "chr22_import": "/home/zmvanw01/projects/EF/240505/run_hifiasm/101-1027/merged_bam/annotations/101-1027/chr22/101-1027_make_gene_file_imported.csv",
        "igh_gene": "/home/zmvanw01/projects/EF/240505/run_hifiasm/101-1027/merged_bam/annotations/101-1027/igh/101-1027_make_gene_file.csv",
        "igh_import": "/home/zmvanw01/projects/EF/240505/run_hifiasm/101-1027/merged_bam/annotations/101-1027/igh/101-1027_make_gene_file_imported.csv"
    }

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split()
            sample = fields[0]
            ighc_gene = fields[1]
            ighc_import = fields[2]
            asm_bam = fields[3]
            ccs_bam = fields[4]

            new_entry = f"{sample}\t{asm_bam}\t{fake_paths['chr2_gene']}\t{fake_paths['chr2_import']}\t{fake_paths['chr22_gene']}\t{fake_paths['chr22_import']}\t{fake_paths['igh_gene']}\t{fake_paths['igh_import']}\t{ighc_gene}\t{ighc_import}\t{ccs_bam}"
            outfile.write(new_entry + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_fofn.py <input_fofn> <output_fofn>")
        sys.exit(1)

    input_fofn = sys.argv[1]
    output_fofn = sys.argv[2]
    process_fofn(input_fofn, output_fofn)
