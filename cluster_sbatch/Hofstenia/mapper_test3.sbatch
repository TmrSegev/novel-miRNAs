#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name mapper			### name of the job
#SBATCH --output mapper-%J.out			### output log for running job - %J for job number
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

echo start ## execute python script – replace with your own command 
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/IST1.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_IST1.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_IST1.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/IST2.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_IST2.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_IST2.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/IST3.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_IST3.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_IST3.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PDi1.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PDi1.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PDi1.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PDi2.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PDi2.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PDi2.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PDi3.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PDi3.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PDi3.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PDii1.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PDii1.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PDii1.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PDii2.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PDii2.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PDii2.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PDii3.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PDii3.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PDii3.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PH1.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PH1.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PH1.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PH2.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PH2.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PH2.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PH3.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PH3.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PH3.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PL1.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PL1.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PL1.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PL2.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PL2.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PL2.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/PL3.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_PL3.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_PL3.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/SMA1.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_SMA1.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_SMA1.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/SMA2.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_SMA2.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_SMA2.fasta
mapper.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/SMA3.filtered.fastq -e -i -j -m -h -p /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/hofsteniaGenomeIndexed -t ../mapper_out_test/hofstenia_Seq_vs_genome_SMA3.arf -s ../mapper_out_test/hofstenia_Seq_collapsed_SMA3.fasta
echo finished.