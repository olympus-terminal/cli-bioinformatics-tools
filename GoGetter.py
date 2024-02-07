import pandas as pd
import sys

from Bio import Entrez
Entrez.email = "drn2@nyu.edu"  # Always tell NCBI who you are

def get_go_terms(pfam_id):
    # Fetching the information about this PFAM from NCBI Entrez, using the Biopython 
    handle = Entrez.esearch(db="protein", term=pfam_id)
    record = Entrez.read(handle)
    handle.close()
    
    # Now we fetch detailed information for each ID we have found
    gos = []
    for i in record["IdList"]:
        handle = Entrez.efetch(db="protein", id=i, rettype="gb", retmode="text")
        for line in handle:
            if 'GO:' in line:
                info = line[line.find('GO:')+3:].strip()
                info_split = info.split(' ')[0]
                gos.append(info_split)
        handle.close()
    return ';'.join(list(set(gos)))

def add_go_terms_infile(input_file, output_file):
    # Load the .csv file.
    df = pd.read_csv(input_file)

    # Add a new column storing the GO terms associated with the PFAM from the 'Accession' column.
    df['GO_terms'] = df['Accession'].apply(get_go_terms)

    # Save the dataframe back to a .csv file.
    df.to_csv(output_file, index=False)

infile = sys.argv[1]
outfile = sys.argv[2]
add_go_terms_infile(infile, outfile)
