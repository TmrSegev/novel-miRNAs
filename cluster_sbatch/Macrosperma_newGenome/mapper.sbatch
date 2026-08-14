#!/bin/bash

#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mapper-macrosperma-v2
#SBATCH --output mapper-%J.out
#SBATCH --mail-user=tzahben@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

# Per-library mapper.pl (Hofstenia-style): one FASTQ per command; no config.txt / -d.
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq
IDX=../genome/index/macrospermaNewGenomeIndexed
OUT=../mapper_out

echo start
mapper.pl ${FQ}/SRR13072564.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR4.arf -s ${OUT}/macrosperma_Seq_collapsed_MR4.fasta
mapper.pl ${FQ}/SRR13072565.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR5.arf -s ${OUT}/macrosperma_Seq_collapsed_MR5.fasta
mapper.pl ${FQ}/SRR13072566.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR6.arf -s ${OUT}/macrosperma_Seq_collapsed_MR6.fasta
mapper.pl ${FQ}/SRR13072567.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR7.arf -s ${OUT}/macrosperma_Seq_collapsed_MR7.fasta
mapper.pl ${FQ}/SRR13072568.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR8.arf -s ${OUT}/macrosperma_Seq_collapsed_MR8.fasta
echo finished.
