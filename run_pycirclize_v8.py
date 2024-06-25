import argparse
from pycirclize import Circos
import pandas as pd
import numpy as np
import re

def parse_arguments():
    parser = argparse.ArgumentParser(description="Visualize genome data using Pycirclize")
    parser.add_argument("karyotype", help="Path to the karyotype file")
    parser.add_argument("coverage", help="Path to the coverage BED file (chr, start, end, coverage)")
    parser.add_argument("genes", help="Path to the TransDecoder gene annotation BED file")
    parser.add_argument("output", help="Path to the output SVG file")
    return parser.parse_args()

def load_karyotype(file_path):
    karyotype = pd.read_csv(file_path, sep="\s+", header=None, 
                            names=["chr", "dash", "id", "id2", "start", "end", "id3"])
    karyotype['size'] = karyotype['end'] - karyotype['start']
    top_25 = karyotype.nlargest(25, 'size')
    return {row["id"]: row["size"] for _, row in top_25.iterrows()}

def calculate_space_size(sectors):
    total_genome_size = sum(sectors.values())
    return max(0.5, min(2, total_genome_size * 0.0000001))

def parse_transdecoder_bed(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 6:
                chr_name = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                info = fields[3]
                
                match = re.search(r'GENE\.([^~]+)', info)
                gene_id = match.group(1) if match else 'Unknown'
                
                data.append({
                    'chr': chr_name,
                    'start': start,
                    'end': end,
                    'gene': gene_id
                })
    
    return pd.DataFrame(data)

def read_coverage_file(file_path):
    def parse_coverage(x):
        return float(x) if x != '.' else np.nan

    coverage = pd.read_csv(file_path, sep="\t", header=None, 
                           names=["chr", "start", "end", "coverage"],
                           dtype={'chr': str, 'start': int, 'end': int, 'coverage': str})
    coverage['coverage'] = coverage['coverage'].apply(parse_coverage)
    return coverage

def main():
    try:
        args = parse_arguments()

        sectors = load_karyotype(args.karyotype)
        space_size = calculate_space_size(sectors)
        
        print(f"Total genome size: {sum(sectors.values()):,}")
        print(f"Number of chromosomes: {len(sectors)}")
        print(f"Calculated space size: {space_size}")

        circos = Circos(sectors, space=space_size)

        coverage = read_coverage_file(args.coverage)
        genes = parse_transdecoder_bed(args.genes)

        for sector in circos.sectors:
            chromosome = sector.name
            
            sector.text(f"Chr {chromosome}", r=110, size=8)
            
            chr_coverage = coverage[coverage["chr"] == chromosome]
            if not chr_coverage.empty:
                coverage_track = sector.add_track((75, 100), r_pad_ratio=0.1)
                coverage_track.axis(fc="#EDEDED")
                coverage_track.xticks_by_interval(int(sector.size / 5))
                
                # Remove NaN values before plotting
                valid_coverage = chr_coverage.dropna(subset=['coverage'])
                if not valid_coverage.empty:
                    coverage_track.line(valid_coverage["start"].values, valid_coverage["coverage"].values, color="blue")
            
            chr_genes = genes[genes["chr"] == chromosome]
            if not chr_genes.empty:
                gene_track = sector.add_track((45, 70), r_pad_ratio=0.1)
                gene_track.axis(fc="#EDEDED")
                for _, gene in chr_genes.iterrows():
                    gene_track.rect((gene["start"], gene["end"]), 1, fc="red")
                    short_gene_id = gene["gene"].split('-')[0]
                    gene_track.text(short_gene_id, gene["start"], r=60, size=6, orientation="vertical")

        circos.savefig(args.output, dpi=300)
        print(f"Visualization saved to {args.output}")

    except ValueError as ve:
        print(f"Error: {ve}")
        print("Please check your input files and ensure they contain valid data.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
