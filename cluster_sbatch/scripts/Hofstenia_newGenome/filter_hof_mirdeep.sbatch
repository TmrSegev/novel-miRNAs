#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name filter			### name of the job
#SBATCH --output filter-hof-mirdeep-%J.out			### output log for running job - %J for job number
##SBATCH --mail-user=stome@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=4G				### ammount of RAM memory
#SBATCH --cpus-per-task=12			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/my_env				### activate a conda environment, replace my_env with your conda environmen

echo start ## execute python script � replace with your own command 
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia_newGenome/mirdeep_out/
for folder in EC1 EC2 EC3 GA1 GA2 GA3 DI1 DI2 DI3 PDi1 PDi2 PDi3 PDii1 PDii2 PDii3 PL1 PL2 PL3 PH1 PH2 PH3 HL1 HL2 HL3 IST1 IST2 IST3 AMP1 AMP2 AMP3 SMA1 SMA2 SMA3; do
    cd "$folder"
    echo "$folder"
    python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepPerLibraryFilter.py -i result*.csv --filter-s 10 --filter-mc 10 --exclude-c 100
    cd ..
done
echo finished.