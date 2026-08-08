#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name mapper_job			### name of the job
#SBATCH --output mapper-%J.out			### output log for running job - %J for job number
#SBATCH --mail-user=stome@post.bgu.ac.il	### user's email for sending job status messages
#SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=128G				### ammount of RAM memory
#SBATCH --cpus-per-task=12			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /sise/home/stome/.conda/envs/my_env/				### activate a conda environment, replace my_env with your conda environmen

# Per-library mapper.pl (Hofstenia-style): one FASTQ per command; no config.txt / -d.
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/TrimmedFastq
IDX=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/Genome/Index/elegansGenomeIndexed
OUT=../mapper_out

echo start
mapper.pl ${FQ}/SRR13072557.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE57.arf -s ${OUT}/elegans_Seq_collapsed_CE57.fasta
mapper.pl ${FQ}/SRR13072558.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE58.arf -s ${OUT}/elegans_Seq_collapsed_CE58.fasta
mapper.pl ${FQ}/SRR13072559.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE59.arf -s ${OUT}/elegans_Seq_collapsed_CE59.fasta
mapper.pl ${FQ}/SRR13072560.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE60.arf -s ${OUT}/elegans_Seq_collapsed_CE60.fasta
mapper.pl ${FQ}/SRR13072561.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE61.arf -s ${OUT}/elegans_Seq_collapsed_CE61.fasta
mapper.pl ${FQ}/SRR13072562.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE62.arf -s ${OUT}/elegans_Seq_collapsed_CE62.fasta
mapper.pl ${FQ}/SRR13072563.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE63.arf -s ${OUT}/elegans_Seq_collapsed_CE63.fasta
mapper.pl ${FQ}/SRR13072569.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE69.arf -s ${OUT}/elegans_Seq_collapsed_CE69.fasta
mapper.pl ${FQ}/SRR13072578.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE78.arf -s ${OUT}/elegans_Seq_collapsed_CE78.fasta
mapper.pl ${FQ}/SRR13072579.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE79.arf -s ${OUT}/elegans_Seq_collapsed_CE79.fasta
mapper.pl ${FQ}/SRR13072580.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE80.arf -s ${OUT}/elegans_Seq_collapsed_CE80.fasta
mapper.pl ${FQ}/SRR13072581.1_trimmed.fastq -e -i -j -m -h -p ${IDX} -t ${OUT}/elegans_Seq_vs_genome_CE81.arf -s ${OUT}/elegans_Seq_collapsed_CE81.fasta
echo finished.
