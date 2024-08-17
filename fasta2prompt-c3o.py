import sys

def process_fasta(fasta_file):
  """
  Processes a FASTA file to remove headers, non-amino acid characters,
  text wrapping, and appends '<' to the end of each sequence.

  Each sequence is processed individually and printed on a separate line.

  Args:
    fasta_file: Path to the FASTA file.
  """

  with open(fasta_file, 'r') as f:
    sequence = ''
    for line in f:
      if line.startswith('>'):  # New sequence
        if sequence:  # Print previous sequence if it exists
          print(sequence + '<') 
        sequence = ''  # Reset for the new sequence
      else:
        # Remove non-amino acid characters and whitespace
        cleaned_line = ''.join(c for c in line.strip() if c.isalpha())
        sequence += cleaned_line

    # Print the last sequence
    if sequence:
      print(sequence + '<')

if __name__ == "__main__":
  if len(sys.argv) != 2:
    print("Usage: python process_fasta.py <fasta_file>")
    sys.exit(1)

  fasta_file = sys.argv[1]
  process_fasta(fasta_file)
