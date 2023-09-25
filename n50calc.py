def calculate_n50(fasta_file):
    """
    Calculate the N50 value of a genomic assembly in FASTA format.
    
    Parameters:
    - fasta_file: Path to the FASTA file.

    Returns:
    - N50 value.
    """
    lengths = []
    total_length = 0

    with open(fasta_file, 'r') as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    lengths.append(len(seq))
                    total_length += len(seq)
                    seq = ""
            else:
                seq += line.strip()
        if seq:  # Add the last sequence
            lengths.append(len(seq))
            total_length += len(seq)

    # Sort lengths in descending order
    sorted_lengths = sorted(lengths, reverse=True)

    # Calculate N50
    cum_length = 0
    for length in sorted_lengths:
        cum_length += length
        if cum_length >= total_length / 2:
            return length

def main(input_file, output_file):
    n50 = calculate_n50(input_file)
    with open(output_file, 'w') as out:
        out.write("Assembly\tN50\n")
        out.write(f"{input_file}\t{n50}\n")

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python n50_calculator_single.py <input_fasta_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
