#!/bin/bash

#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mapper-elegans-wbps19
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
# Reads reused from Elegans/TrimmedFastq; index is WBPS19.
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/TrimmedFastq
IDX=../genome/index/elegansNewGenomeIndexed
OUT=../mapper_out

echo start
mapper.pl ${FQ}/SRR13072557.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE57.arf -s ${OUT}/elegans_Seq_collapsed_CE57.fasta
mapper.pl ${FQ}/SRR13072558.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE58.arf -s ${OUT}/elegans_Seq_collapsed_CE58.fasta
mapper.pl ${FQ}/SRR13072559.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE59.arf -s ${OUT}/elegans_Seq_collapsed_CE59.fasta
mapper.pl ${FQ}/SRR13072560.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE60.arf -s ${OUT}/elegans_Seq_collapsed_CE60.fasta
mapper.pl ${FQ}/SRR13072561.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE61.arf -s ${OUT}/elegans_Seq_collapsed_CE61.fasta
mapper.pl ${FQ}/SRR13072562.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE62.arf -s ${OUT}/elegans_Seq_collapsed_CE62.fasta
mapper.pl ${FQ}/SRR13072563.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE63.arf -s ${OUT}/elegans_Seq_collapsed_CE63.fasta
mapper.pl ${FQ}/SRR13072569.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE69.arf -s ${OUT}/elegans_Seq_collapsed_CE69.fasta
mapper.pl ${FQ}/SRR13072578.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE78.arf -s ${OUT}/elegans_Seq_collapsed_CE78.fasta
mapper.pl ${FQ}/SRR13072579.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE79.arf -s ${OUT}/elegans_Seq_collapsed_CE79.fasta
mapper.pl ${FQ}/SRR13072580.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE80.arf -s ${OUT}/elegans_Seq_collapsed_CE80.fasta
mapper.pl ${FQ}/SRR13072581.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE81.arf -s ${OUT}/elegans_Seq_collapsed_CE81.fasta
echo finished.
