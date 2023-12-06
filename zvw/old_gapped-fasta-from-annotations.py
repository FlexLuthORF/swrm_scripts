import csv
import sys

def process_csv_to_fasta(csv_file):
    fasta_data = {'V': [], 'V_gapped': [], 'D': [], 'J': []}
    fasta_filenames = {
        'V': 'human_ig_V_alleles.fasta',
        'V_gapped': 'human_ig_V_alleles_gapped.fasta',
        'D': 'human_ig_D_alleles.fasta',
        'J': 'human_ig_J_alleles.fasta'
    }

    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            allele = row['vdjbase_allele']
            if len(allele) >= 4:
                category = allele[3]
                if category in ['V', 'D', 'J']:
                    region = row.get(f'{category}-REGION', '')
                    if region:
                        fasta_data[category].append(f'>{allele}\n{region}')

                    if category == 'V':  # Additionally process V-REGION-GAPPED
                        region_gapped = row.get('V-REGION-GAPPED', '')
                        if region_gapped:
                            fasta_data['V_gapped'].append(f'>{allele}\n{region_gapped}')
            else:
                print(f"Skipped allele: {allele}")  # Debug print

    for key, data in fasta_data.items():
        if data:
            with open(fasta_filenames[key], 'w') as fasta_file:
                fasta_file.write('\n'.join(data))
        else:
            print(f"No data for category {key}, no file created.")  # Debug print

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <csv_file>")
    else:
        process_csv_to_fasta(sys.argv[1])
