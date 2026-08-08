#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name my_job			### name of the job
#SBATCH --output mapper-%J.out			### output log for running job - %J for job number
#SBATCH --mail-user=tzahben@post.bgu.ac.il	### user's email for sending job status messages
#SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=128G				### ammount of RAM memory
#SBATCH --cpus-per-task=12			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/my_env				### activate a conda environment, replace my_env with your conda environmen

# Per-library mapper.pl (Hofstenia-style): one FASTQ per command; no config.txt / -d.
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq
IDX=../genome/index/macrospermaGenomeIndexed
OUT=../mapper_out

echo start
mapper.pl ${FQ}/SRR13072564.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR4.arf -s ${OUT}/macrosperma_Seq_collapsed_MR4.fasta
mapper.pl ${FQ}/SRR13072565.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR5.arf -s ${OUT}/macrosperma_Seq_collapsed_MR5.fasta
mapper.pl ${FQ}/SRR13072566.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR6.arf -s ${OUT}/macrosperma_Seq_collapsed_MR6.fasta
mapper.pl ${FQ}/SRR13072567.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR7.arf -s ${OUT}/macrosperma_Seq_collapsed_MR7.fasta
mapper.pl ${FQ}/SRR13072568.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/macrosperma_Seq_vs_genome_MR8.arf -s ${OUT}/macrosperma_Seq_collapsed_MR8.fasta
echo finished.
