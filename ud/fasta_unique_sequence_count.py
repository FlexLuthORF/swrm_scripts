import sys

from Bio import SeqIO
from collections import Counter

def count_unique_sequences(fasta_file):
    sequences = []

    # Iterate through each sequence in the FASTA file and collect sequences
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))

    # Count occurrences of each sequence
    sequence_counter = Counter(sequences)

    # Sort sequences by count in descending order and print top 4
    #print("Top 4 most frequent sequences:")
    for sequence, count in sequence_counter.most_common(5):
        #print(f"{sequence}\t{count}")
        print(f"{count}")
if __name__ == "__main__":
    fasta_file = sys.argv[1]  # Specify input FASTA file
    count_unique_sequences(fasta_file)
