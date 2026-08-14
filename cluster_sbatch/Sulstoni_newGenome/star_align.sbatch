#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name STAR-align-Sulstoni_newGenome
#SBATCH --output STAR-align-%J.out
#SBATCH --mail-user=tzahben@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda

source activate /home/stome/.conda/envs/my_env

TRIM=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/TrimmedFastq
echo start
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072570.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR0/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072571.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR1/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072572.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR2/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072573.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR3/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072574.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR4/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072575.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR5/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072576.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR6/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072577.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/SR7/Sulstoni_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
echo finished.
