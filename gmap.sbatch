#!/bin/bash

#SBATCH --mem=90GB
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=28

gmap_build -d genome -D ./ ./ctg.fa

wait

gmap -d genome -D ./ -f 2 -n 0 -t 28 \
    transdecoder.cds > gmap_aligned_cds.gff3
