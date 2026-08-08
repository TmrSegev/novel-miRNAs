#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name my_job			### name of the job
#SBATCH --output blast_macrosperma-%J.out			### output log for running job - %J for job number
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
blastn -query ../../Charles_seq/Macrosperma/scripts/Macrosperma_mirdeep.fasta -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB -out ../queries/Macrosperma/miRdeep_blastn_compact -outfmt 6 -evalue 10 -task blastn-short
blastn -query /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts/Macrosperma_sRNAbench.fasta -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB -out ../queries/Macrosperma/sRNAbench_blastn_compact -outfmt 6 -evalue 10 -task blastn-short
echo finished.


