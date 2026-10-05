import random
import argparse

def introduce_errors(sequence, error_rate):
    """
    Introduce errors into a sequence based on the given error rate.

    :param sequence: Original sequence
    :param error_rate: Probability of introducing an error at each position
    :return: Sequence with errors introduced
    """
    nucleotides = ['A', 'C', 'G', 'T']
    mutated_sequence = list(sequence)

    for i in range(len(sequence)):
        if random.random() < error_rate:
            original_nucleotide = sequence[i]
            possible_mutations = [nuc for nuc in nucleotides if nuc != original_nucleotide]
            mutated_sequence[i] = random.choice(possible_mutations)

    return ''.join(mutated_sequence)

def introduce_mixed_errors(sequence, error_rate):
    """
    Introduce errors as substitutions, insertions and deletions in equal proportion.

    :param sequence: Original sequence
    :param error_rate: Probability of introducing an error at each position
    :return: Sequence with errors introduced
    """
    nucleotides = ['A', 'C', 'G', 'T']
    mutated_sequence = []

    for nucleotide in sequence:
        if random.random() >= error_rate:
            mutated_sequence.append(nucleotide)
            continue
        kind = random.randrange(3)
        if kind == 0:  # substitution
            mutated_sequence.append(random.choice([nuc for nuc in nucleotides if nuc != nucleotide]))
        elif kind == 1:  # insertion
            mutated_sequence.append(nucleotide)
            mutated_sequence.append(random.choice(nucleotides))
        # otherwise a deletion

    return ''.join(mutated_sequence)

def generate_telomere_sequences(motif, error_rate, sequence_length, num_sequences, error_model="substitution"):
    """
    Generate telomere sequences with specified parameters.

    :param motif: Telomere motif
    :param error_rate: Error rate for mutations
    :param sequence_length: Length of each generated sequence
    :param num_sequences: Number of sequences to generate
    :param error_model: "substitution" or "mixed" (substitutions and indels)
    :return: List of generated telomere sequences
    """
    motif_length = len(motif)
    sequences = []

    for _ in range(num_sequences):
        sequence = ''
        while len(sequence) < sequence_length:
            sequence += motif
        sequence = sequence[:sequence_length]  # Trim to desired length

        if error_model == "mixed":
            mutated_sequence = introduce_mixed_errors(sequence, error_rate)
        else:
            mutated_sequence = introduce_errors(sequence, error_rate)
        sequences.append(mutated_sequence)

    return sequences

def save_to_fasta(sequences, error_rate, sequence_length, output_file):
    """
    Save sequences to a FASTA file with headers.

    :param sequences: List of sequences to save
    :param error_rate: Error rate used in generating sequences
    :param sequence_length: Length of each sequence
    :param output_file: Path to the output file
    """
    with open(output_file, 'w') as f:
        for i, seq in enumerate(sequences):
            header = f">seq{i+1}_len{sequence_length}_err{error_rate}\n"
            f.write(header)
            f.write(seq + '\n')

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate telomere sequences with specified error rate.')
    parser.add_argument('motif', type=str, help='Telomere motif')
    parser.add_argument('error_rate', type=float, help='Error rate (e.g., 0.01 for 1 percent)')
    parser.add_argument('sequence_length', type=int, help='Length of each generated sequence')
    parser.add_argument('num_sequences', type=int, help='Number of sequences to generate')
    parser.add_argument('output_file', type=str, help='Output file to save sequences in FASTA format')
    parser.add_argument('--error-model', choices=['substitution', 'mixed'], default='substitution',
                        help='Substitutions only (as in the paper), or substitutions and indels in equal proportion')

    args = parser.parse_args()

    # Generate sequences
    sequences = generate_telomere_sequences(args.motif, args.error_rate, args.sequence_length, args.num_sequences, args.error_model)

    # Save sequences to FASTA file
    save_to_fasta(sequences, args.error_rate, args.sequence_length, args.output_file)

if __name__ == "__main__":
    main()
