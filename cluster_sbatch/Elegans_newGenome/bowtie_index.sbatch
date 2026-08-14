#!/bin/bash

#SBATCH --partition main
#SBATCH --time 1-00:00:00
#SBATCH --job-name bowtie-elegans-wbps19
#SBATCH --output bowtie_index-%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

GENOME="../genome/CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna"
bowtie-build -f "$GENOME" ../genome/index/elegansNewGenomeIndexed
samtools faidx "$GENOME"

echo finished.
