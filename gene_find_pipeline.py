#!/bin/bash

# HPC-ready HISAT2 + StringTie + BRAKER Pipeline for RNA-seq Mapping and Gene Finding

# Function to display usage information
usage() {
    echo "Usage: $0 -g <reference_genome> -1 <forward_reads> -2 <reverse_reads> -i <ab_initio_gff3> -o <output_prefix> [-t <threads>]"
    echo "  -g: Reference genome file (FASTA)"
    echo "  -1: Forward reads file (FASTQ)"
    echo "  -2: Reverse reads file (FASTQ)"
    echo "  -i: Ab initio gene predictions (GFF3)"
    echo "  -o: Output prefix for all generated files"
    echo "  -t: Number of threads to use (default: 8)"
    exit 1
}

# Parse command line arguments
while getopts "g:1:2:i:o:t:" opt; do
    case $opt in
        g) GENOME=$OPTARG ;;
        1) FWD_READS=$OPTARG ;;
        2) REV_READS=$OPTARG ;;
        i) AB_INITIO_GFF3=$OPTARG ;;
        o) OUTPUT_PREFIX=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        *) usage ;;
    esac
done

# Check if required arguments are provided
if [ -z "$GENOME" ] || [ -z "$FWD_READS" ] || [ -z "$REV_READS" ] || [ -z "$AB_INITIO_GFF3" ] || [ -z "$OUTPUT_PREFIX" ]; then
    usage
fi

# Set default number of threads if not provided
THREADS=${THREADS:-8}

# Step 1: Index the reference genome with HISAT2
echo "Indexing reference genome..."
hisat2-build $GENOME ${OUTPUT_PREFIX}_genome_index

# Step 2: Align RNA-seq reads to the reference genome
echo "Aligning reads to reference genome..."
hisat2 -p $THREADS --dta -x ${OUTPUT_PREFIX}_genome_index -1 $FWD_READS -2 $REV_READS -S ${OUTPUT_PREFIX}.sam

# Step 3: Convert SAM to BAM, sort, and index
echo "Converting SAM to BAM, sorting, and indexing..."
samtools view -bS ${OUTPUT_PREFIX}.sam > ${OUTPUT_PREFIX}.bam
samtools sort ${OUTPUT_PREFIX}.bam -o ${OUTPUT_PREFIX}.sorted.bam
samtools index ${OUTPUT_PREFIX}.sorted.bam

# Step 4: Assemble transcripts with StringTie
echo "Assembling transcripts..."
stringtie ${OUTPUT_PREFIX}.sorted.bam -o ${OUTPUT_PREFIX}.gtf -p $THREADS

# Step 5: Convert StringTie GTF to GFF3
echo "Converting StringTie GTF to GFF3..."
gffread ${OUTPUT_PREFIX}.gtf -o ${OUTPUT_PREFIX}_stringtie.gff3

# Step 6: Merge ab initio predictions with RNA-seq evidence
echo "Merging ab initio predictions with RNA-seq evidence..."
bedtools intersect -a $AB_INITIO_GFF3 -b ${OUTPUT_PREFIX}_stringtie.gff3 -wa -u > ${OUTPUT_PREFIX}_merged_evidence.gff3

# Step 7: Run BRAKER for gene prediction
echo "Running BRAKER for gene prediction..."
braker.pl --genome=$GENOME --bam=${OUTPUT_PREFIX}.sorted.bam --gff3=${OUTPUT_PREFIX}_merged_evidence.gff3 --cores=$THREADS --workingdir=${OUTPUT_PREFIX}_braker_out

# Step 8: Extract final gene predictions
echo "Extracting final gene predictions..."
cp ${OUTPUT_PREFIX}_braker_out/augustus.hints.gff3 ${OUTPUT_PREFIX}_final_gene_predictions.gff3

echo "Pipeline completed. Final gene predictions are in '${OUTPUT_PREFIX}_final_gene_predictions.gff3'."
