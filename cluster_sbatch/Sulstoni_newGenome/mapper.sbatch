#!/bin/bash

#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mapper-sulstoni-wbps19
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
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/TrimmedFastq
IDX=../genome/index/sulstoniNewGenomeIndexed
OUT=../mapper_out

echo start
mapper.pl ${FQ}/SRR13072570.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR0.arf -s ${OUT}/sulstoni_Seq_collapsed_SR0.fasta
mapper.pl ${FQ}/SRR13072571.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR1.arf -s ${OUT}/sulstoni_Seq_collapsed_SR1.fasta
mapper.pl ${FQ}/SRR13072572.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR2.arf -s ${OUT}/sulstoni_Seq_collapsed_SR2.fasta
mapper.pl ${FQ}/SRR13072573.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR3.arf -s ${OUT}/sulstoni_Seq_collapsed_SR3.fasta
mapper.pl ${FQ}/SRR13072574.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR4.arf -s ${OUT}/sulstoni_Seq_collapsed_SR4.fasta
mapper.pl ${FQ}/SRR13072575.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR5.arf -s ${OUT}/sulstoni_Seq_collapsed_SR5.fasta
mapper.pl ${FQ}/SRR13072576.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR6.arf -s ${OUT}/sulstoni_Seq_collapsed_SR6.fasta
mapper.pl ${FQ}/SRR13072577.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR7.arf -s ${OUT}/sulstoni_Seq_collapsed_SR7.fasta
echo finished.
