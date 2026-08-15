#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name STAR-align-Elegans_newGenome
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

TRIM=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/TrimmedFastq
mkdir -p ../STAR/align_to_genome/{CE57,CE58,CE59,CE60,CE61,CE62,CE63,CE69,CE78,CE79,CE80,CE81}
echo start
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072557.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE57/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072558.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE58/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072559.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE59/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072560.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE60/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072561.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE61/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072562.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE62/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072563.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE63/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072569.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE69/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072578.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE78/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072579.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE79/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072580.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE80/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn ${TRIM}/SRR13072581.1_trimmed.fastq --outFileNamePrefix ../STAR/align_to_genome/CE81/Elegans_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
echo finished.
