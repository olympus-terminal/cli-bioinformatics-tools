import csv
import sys
from pathlib import Path

def csv_to_fasta(input_file, output_file):
    """
    Convert a CSV file containing protein sequences to FASTA format.
    
    Args:
        input_file (str): Path to input CSV file
        output_file (str): Path to output FASTA file
    """
    try:
        with open(input_file, 'r') as csv_file, open(output_file, 'w') as fasta_file:
            # Use csv reader with flexible handling of whitespace
            csv_reader = csv.reader(csv_file, skipinitialspace=True)
            
            for row in csv_reader:
                if len(row) < 2:
                    print(f"Warning: Skipping malformed line in {input_file}")
                    continue
                    
                # Get identifier and sequence
                identifier = row[0].strip()
                sequence = row[1].strip()
                
                # Write in FASTA format
                fasta_file.write(f">{identifier}\n")
                
                # Write sequence in lines of 60 characters
                for i in range(0, len(sequence), 60):
                    fasta_file.write(f"{sequence[i:i+60]}\n")
                    
        print(f"Successfully converted {input_file} to {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Could not find input file {input_file}")
        return False
    except PermissionError:
        print(f"Error: Permission denied when accessing {input_file} or {output_file}")
        return False
    except Exception as e:
        print(f"Error: An unexpected error occurred: {str(e)}")
        return False
    
    return True

def main():
    # Process all CSV files with 'noAst.csv' suffix in current directory
    csv_files = Path('.').glob('*noAst.csv')
    
    for csv_file in csv_files:
        # Create output filename by replacing .csv extension with .fasta
        output_file = csv_file.stem.replace('.csv', '') + '.fasta'
        
        # Convert file
        if csv_to_fasta(str(csv_file), output_file):
            print(f"Processed {csv_file} -> {output_file}")
        else:
            print(f"Failed to process {csv_file}")

if __name__ == "__main__":
    main()
