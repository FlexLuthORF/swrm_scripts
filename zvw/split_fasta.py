import sys

def split_fasta_by_loci(input_fasta):
    # Derive output file names from input file name
    base_name = input_fasta.rsplit('.', 1)[0]
    j_output_filename = f'{base_name}_J_loci.fasta'
    v_output_filename = f'{base_name}_V_loci.fasta'
    
    # Initialize file handles outside the with statement to avoid UnboundLocalError
    j_file = None
    v_file = None

    try:
        # Open the input file for reading
        with open(input_fasta, 'r') as fasta:
            for line in fasta:
                # If it's a header line, decide which file to write to based on the fourth character (excluding '>')
                if line.startswith('>'):
                    if line[4] == 'J':
                        if not j_file:
                            j_file = open(j_output_filename, 'w')
                        current_file = j_file
                    elif line[4] == 'V':
                        if not v_file:
                            v_file = open(v_output_filename, 'w')
                        current_file = v_file
                    else:
                        current_file = None
                # Write the line (either a header or a sequence line) to the appropriate file
                if current_file:
                    current_file.write(line)
    finally:
        # Close the output files if they were opened
        if j_file:
            j_file.close()
        if v_file:
            v_file.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py input_fasta.fasta")
    else:
        input_fasta = sys.argv[1]
        split_fasta_by_loci(input_fasta)
