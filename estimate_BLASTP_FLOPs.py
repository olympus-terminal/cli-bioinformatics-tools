def estimate_blastp_flops(query_length, database_size, avg_sequence_length):
    """
    Estimate FLOPs for BLASTP algorithm.
    
    :param query_length: Length of the query sequence
    :param database_size: Number of sequences in the database
    :param avg_sequence_length: Average length of sequences in the database
    :return: Estimated number of FLOPs
    """
    # Constants (these are approximations and may vary based on implementation)
    WORD_SIZE = 3  # Typical word size for proteins
    SCORING_FLOPS = 10  # Approximate FLOPs for scoring a single residue pair
    EXTENSION_FLOPS = 50  # Approximate FLOPs for extending a hit
    
    # Step 1: Word matching
    word_matches = query_length * database_size * avg_sequence_length / (20**WORD_SIZE)
    word_matching_flops = word_matches * WORD_SIZE * SCORING_FLOPS
    
    # Step 2: Ungapped extension
    ungapped_extensions = word_matches * 0.1  # Assume 10% of word matches are extended
    ungapped_extension_flops = ungapped_extensions * 20 * EXTENSION_FLOPS
    
    # Step 3: Gapped extension
    gapped_extensions = ungapped_extensions * 0.1  # Assume 10% of ungapped extensions are further extended
    gapped_extension_flops = gapped_extensions * avg_sequence_length * EXTENSION_FLOPS
    
    total_flops = word_matching_flops + ungapped_extension_flops + gapped_extension_flops
    
    return int(total_flops)

def main():
    # Example usage
    query_length = 300
    database_size = 1000000
    avg_sequence_length = 400

    flops = estimate_blastp_flops(query_length, database_size, avg_sequence_length)
    print(f"Estimated FLOPs for BLASTP: {flops:,}")

    # Interactive mode
    while True:
        try:
            query_length = int(input("Enter query sequence length (or 0 to exit): "))
            if query_length == 0:
                break
            database_size = int(input("Enter database size (number of sequences): "))
            avg_sequence_length = int(input("Enter average sequence length in database: "))

            flops = estimate_blastp_flops(query_length, database_size, avg_sequence_length)
            print(f"Estimated FLOPs for BLASTP: {flops:,}")
        except ValueError:
            print("Please enter valid integer values.")
        print()  # Empty line for readability

if __name__ == "__main__":
    main()
