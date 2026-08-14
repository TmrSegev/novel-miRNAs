#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 6-23:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name STAR3		### name of the job
#SBATCH --output STAR-align3-%J.out			### output log for running job - %J for job number
##SBATCH --mail-user=stome@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=128G				### ammount of RAM memory
#SBATCH --cpus-per-task=12			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/my_env				### activate a conda environment, replace my_env with your conda environmen

echo start ## execute python script – replace with your own command 
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PH1.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/PH1/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PH2.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/PH2/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PH3.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/PH3/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PL1.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/PL1/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PL2.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/PL2/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PL3.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/PL3/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/SMA1.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/SMA1/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/SMA2.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/SMA2/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
STAR --genomeDir ../STAR/genome_index/ --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/SMA3.filtered.fastq --outFileNamePrefix ../STAR/align_to_genome/SMA3/Hofstenia_ --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
echo finished.