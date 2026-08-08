#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name mirge			### name of the job
#SBATCH --output mirge-%J.out			### output log for running job - %J for job number
##SBATCH --mail-user=stome@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=4G				### ammount of RAM memory
#SBATCH --cpus-per-task=8			### number of CPU cores

echo "Job start"
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/mirge_env				### activate a conda environment, replace my_env with your conda environmen
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/

# Loop through all files ending in *.filtered.fastq in the specified folder
#for sample in /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/GA3.filtered.fastq
for sample in /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/*.filtered.fastq
do
  echo "Processing $sample"
  miRge3.0 \
    -s $sample \
    -lib /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/ \
    -on Hofstenia \
    -db miRBase \
    -o /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/results \
    -gff \
    -bam \
    -ie \
    -cpu 8
done


echo "Job finished"