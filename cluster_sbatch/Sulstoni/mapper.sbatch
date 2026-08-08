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
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/TrimmedFastq
IDX=../genome/index/sulstoniGenomeIndexed
OUT=../mapper_out

echo start
mapper.pl ${FQ}/SRR13072570.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR0.arf -s ${OUT}/sulstoni_Seq_collapsed_SR0.fasta
mapper.pl ${FQ}/SRR13072571.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR1.arf -s ${OUT}/sulstoni_Seq_collapsed_SR1.fasta
mapper.pl ${FQ}/SRR13072572.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR2.arf -s ${OUT}/sulstoni_Seq_collapsed_SR2.fasta
mapper.pl ${FQ}/SRR13072573.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR3.arf -s ${OUT}/sulstoni_Seq_collapsed_SR3.fasta
mapper.pl ${FQ}/SRR13072574.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR4.arf -s ${OUT}/sulstoni_Seq_collapsed_SR4.fasta
mapper.pl ${FQ}/SRR13072575.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR5.arf -s ${OUT}/sulstoni_Seq_collapsed_SR5.fasta
mapper.pl ${FQ}/SRR13072576.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR6.arf -s ${OUT}/sulstoni_Seq_collapsed_SR6.fasta
mapper.pl ${FQ}/SRR13072577.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/sulstoni_Seq_vs_genome_SR7.arf -s ${OUT}/sulstoni_Seq_collapsed_SR7.fasta
echo finished.
