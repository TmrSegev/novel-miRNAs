#!/bin/bash

#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name STAR-elegans-wbps19
#SBATCH --output STAR-index-%J.out
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

mkdir -p ../STAR/genome_index
STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir ../STAR/genome_index/ \
  --genomeSAindexNbases 11 \
  --genomeFastaFiles ../genome/CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna

echo finished.
