#!/bin/bash

#SBATCH -o slurm-logs/arrayJob_%A_%a.out
#SBATCH -e slurm-logs/arrayJob_%A_%a.err
#SBATCH -a 1
#SBATCH -p bigmem
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000GB
#SBATCH --time=96:00:00

#LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p filelist.txt)
#echo $LINE

#meryl count k=21 GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap1.p_ctg.fa output GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap1.p_ctg.meryl &
#meryl count k=21 GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap2.p_ctg.fa output GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap2.p_ctg.meryl &

#wait

#meryl difference GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap1.p_ctg.meryl GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap2.p_ctg.meryl output hap1.meryl &
#meryl difference GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap2.p_ctg.meryl GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap1.p_ctg.meryl output hap2.meryl 

#meryl count k=21 /scratch/drn2/PROJECTS/GORILLA/Fin-2batchHIFI-plus-ONT/MtnGorillaREADs/Gorilla_1st2ndHIFI.fq output hifi_reads.meryl
#meryl count k=21 /scratch/drn2/PROJECTS/GORILLA/Fin-2batchHIFI-plus-ONT/MtnGorillaREADs/GORILLA-ONT_reads.fq output ONT_reads.meryl

$MERQURY/merqury.sh hifi_reads.meryl hap1.meryl hap2.meryl GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap1.p_ctg.fa GORILLA-ONT_reads_2batch-all-HIFI_ONT.asm.bp.hap1.p_ctg.fa merqury-hifi
