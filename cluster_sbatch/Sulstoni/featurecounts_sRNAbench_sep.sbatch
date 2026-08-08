#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name featurecounts_sRNAbench_sep			### name of the job
#SBATCH --output featurecounts_sRNAbenchp_sep-%J.out			### output log for running job - %J for job number
##SBATCH --mail-user=stome@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=16G				### ammount of RAM memory
#SBATCH --cpus-per-task=12			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/my_env				### activate a conda environment, replace my_env with your conda environmen
echo start ## execute python script – replace with your own command 
featureCounts -t miRNA -g ID -O -s 1 -M -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts/Sulstoni_sRNAbench.gff3 -o ../counts_sep/miRNA_sRNAbench_counts.txt ../STAR/align_to_genome/SR7/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR6/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR5/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR4/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR3/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR2/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR1/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR0/Sulstoni_Aligned.out.sam
echo finished.


