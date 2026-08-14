#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name STAR-align-Macrosperma_newGenome
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

TRIM=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq
echo start
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072564.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/MR4/Macrosperma_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072565.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/MR5/Macrosperma_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072566.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/MR6/Macrosperma_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072567.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/MR7/Macrosperma_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072568.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/MR8/Macrosperma_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
echo finished.
