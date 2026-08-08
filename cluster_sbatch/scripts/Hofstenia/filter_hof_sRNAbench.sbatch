#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name filter			### name of the job
#SBATCH --output filter-hof-bench-%J.out			### output log for running job - %J for job number
##SBATCH --mail-user=stome@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=4G				### ammount of RAM memory
#SBATCH --cpus-per-task=8			### number of CPU cores
### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/my_env				### activate a conda environment, replace my_env with your conda environmen

echo start
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out
for folder in Hofstenia_EC1 Hofstenia_EC2 Hofstenia_EC3 Hofstenia_GA1 Hofstenia_GA2 Hofstenia_GA3 Hofstenia_DI1 Hofstenia_DI2 Hofstenia_DI3 Hofstenia_PDi1 Hofstenia_PDi2 Hofstenia_PDi3 Hofstenia_PDii1 Hofstenia_PDii2 Hofstenia_PDii3 Hofstenia_PL1 Hofstenia_PL2 Hofstenia_PL3 Hofstenia_PH1 Hofstenia_PH2 Hofstenia_PH3 Hofstenia_HL1 Hofstenia_HL2 Hofstenia_HL3 Hofstenia_IST1 Hofstenia_IST2 Hofstenia_IST3 Hofstenia_AMP1 Hofstenia_AMP2 Hofstenia_AMP3 Hofstenia_SMA1 Hofstenia_SMA2 Hofstenia_SMA3; do
    cd "$folder"
    echo "$folder"
    python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py -i novel.txt -a novel451.txt --filter-mc 10
    cd ..
done
echo finished.
