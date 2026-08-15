#!/bin/bash

#SBATCH --partition main
#SBATCH --time 1-00:00:00
#SBATCH --job-name bowtie-hofstenia
#SBATCH --output bowtie_index-%J.out
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

GENOME=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/Hmia.030120.fasta
IDX=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed
mkdir -p "$(dirname "$IDX")"
bowtie-build -f "$GENOME" "$IDX"
samtools faidx "$GENOME"
mkdir -p ../../sRNAtoolboxDB/index
cp -a "${IDX}".* ../../sRNAtoolboxDB/index/

echo finished.
