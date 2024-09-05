import os
from collections import Counter


def read_fasta(file_path):
    """Reads a FASTA file and returns the sequences as a single concatenated string."""
    sequences = []
    with open(file_path, 'r') as fasta_file:
        for line in fasta_file:
            if not line.startswith('>'):  # Skip the header lines
                sequences.append(line.strip())
    return ''.join(sequences)


def count_amino_acids(sequence):
    """Counts the frequency of each amino acid in the given sequence."""
    return Counter(sequence)

def get_frequency(aa_counts: Counter):
    """Counts the frequency of the given amino acid in the sequence."""
    total_aa = sum(aa_counts.values())
    aa_freqs = {aa: count / total_aa for aa, count in aa_counts.items()}
    return aa_freqs


def print_amino_acid_frequencies(frequencies):
    """Prints the amino acid frequencies in a readable format."""
    for aa, count in frequencies.items():
        print(f"'{aa}': {round(count, 7)},")


if __name__ == "__main__":
    # Replace with your FASTA file path

    # read all files in fastas
    for file in os.listdir('fastas'):
        if file.endswith('.fasta'):

            file_path = os.path.join('fastas', file)

            # Read the sequences from the FASTA file
            sequence = read_fasta(file_path)

            # Count amino acid frequencies
            counts = count_amino_acids(sequence)

            # Get the frequency of each amino acid
            frequencies = get_frequency(counts)

            print(f'{file.split(".")[0]} = ' + '{')

            # Print the results
            print_amino_acid_frequencies(frequencies)

            print('}')
