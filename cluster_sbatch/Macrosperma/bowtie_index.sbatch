#!/bin/bash

#SBATCH --partition main
#SBATCH --time 1-00:00:00
#SBATCH --job-name bowtie-macrosperma
#SBATCH --output bowtie_index-%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

GENOME="../genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna"
IDX="../genome/index/macrospermaGenomeIndexed"
mkdir -p ../genome/index
bowtie-build -f "$GENOME" "$IDX"
samtools faidx "$GENOME"
mkdir -p ../../sRNAtoolboxDB/index
cp -a "${IDX}".* ../../sRNAtoolboxDB/index/

echo finished.
