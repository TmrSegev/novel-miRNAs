#!/bin/bash

#SBATCH --partition main
#SBATCH --time 1-00:00:00
#SBATCH --job-name bowtie-elegans
#SBATCH --output bowtie_index-%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

GENOME="../Genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa"
IDX="../Genome/Index/elegansGenomeIndexed"
mkdir -p ../Genome/Index
bowtie-build -f "$GENOME" "$IDX"
samtools faidx "$GENOME"
cp -r ../Genome/Index/. ../../sRNAtoolboxDB/index/

echo finished.
