#!/bin/bash
# LEGACY: sequential all-libraries. Prefer sRNAbench_{LIBRARY}.sbatch (one Slurm job per library).
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name sRNAbench-Macrosperma
#SBATCH --output sRNAbench-%J.out
#SBATCH --mail-user=tzahben@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda

source activate /home/stome/.conda/envs/my_env

DB=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB
TRIM=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq
OUT=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Macrosperma
mkdir -p "$OUT"
echo start
java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072564.1_trimmed.fastq \
  output=${OUT}/Macrosperma_MR4 predict=true species=macrospermaGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072565.1_trimmed.fastq \
  output=${OUT}/Macrosperma_MR5 predict=true species=macrospermaGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072566.1_trimmed.fastq \
  output=${OUT}/Macrosperma_MR6 predict=true species=macrospermaGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072567.1_trimmed.fastq \
  output=${OUT}/Macrosperma_MR7 predict=true species=macrospermaGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072568.1_trimmed.fastq \
  output=${OUT}/Macrosperma_MR8 predict=true species=macrospermaGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

echo finished.
