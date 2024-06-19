import os
import subprocess
import sys

# Input files
assembly_file = sys.argv[1]
reads1_file = sys.argv[2]
reads2_file = sys.argv[3]
genes_file = sys.argv[4]

# Output directory
output_dir = "circos_input"

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# 1. Prepare the karyotype file

# Extract chromosome lengths from the assembly FASTA
chromosome_lengths = {}
with open(assembly_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            chromosome_name = line[1:].strip()
            chromosome_lengths[chromosome_name] = 0
        else:
            chromosome_lengths[chromosome_name] += len(line.strip())

# Write the karyotype file
karyotype_file = os.path.join(output_dir, "karyotype.txt")
with open(karyotype_file, 'w') as f:
    for chromosome_name, length in chromosome_lengths.items():
        f.write(f"chr - {chromosome_name} {chromosome_name} 0 {length} {chromosome_name}\n")

# 2. Calculate and prepare coverage data

def calculate_coverage(sam_file, output_file):
    """Calculates coverage from a SAM file and writes to output_file."""
    coverage = {}
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith("@"):
                continue  # Skip header lines
            fields = line.strip().split("\t")
            chrom = fields[2]
            start = int(fields[3])
            end = start + len(fields[9])  # Calculate read end position
            for i in range(start, end):
                coverage.setdefault(chrom, {}).setdefault(i, 0)
                coverage[chrom][i] += 1
    with open(output_file, 'w') as outfile:
        for chrom, pos_counts in coverage.items():
            for pos, count in pos_counts.items():
                outfile.write(f"{chrom}\t{pos}\t{count}\n")

coverage1_file = os.path.join(output_dir, "coverage1.txt")
coverage2_file = os.path.join(output_dir, "coverage2.txt")

calculate_coverage(reads1_file, coverage1_file)
calculate_coverage(reads2_file, coverage2_file)

# 3. Prepare gene data for Circos

# Assuming the genes.bed file has at least 4 columns: chromosome, start, end, gene_name
genes_circos_file = os.path.join(output_dir, "genes.txt")
with open(genes_file, 'r') as infile, open(genes_circos_file, 'w') as outfile:
    for line in infile:
        chrom, start, end, name, *rest = line.strip().split("\t")
        outfile.write(f"hs1 {chrom} {start} {end} {name}\n")  # Adjust 'hs1' if needed

# 4. Create the Circos configuration file

config_file = os.path.join(output_dir, "circos.conf")
with open(config_file, 'w') as f:
    f.write(f"""


karyotype = {karyotype_file}

<ideogram>
<spacing>
default = 0.005r
</spacing>

# Ideogram appearance
thickness = 20p
fill = yes
stroke_color = black
stroke_thickness = 2p
</ideogram>

<plots>

# Plot for genes
<plot>
type = tile
file = {genes_circos_file}
r0 = 0.85r
r1 = 0.95r
stroke_thickness = 0
</plot>

# Plot for read coverage (example, adapt as needed)
<plot>
type = histogram
file = {reads1_file}.bam
r0 = 0.70r
r1 = 0.80r
fill_color = blue
</plot>

</plots>
""")
