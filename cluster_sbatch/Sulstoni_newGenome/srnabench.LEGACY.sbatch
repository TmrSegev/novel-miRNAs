#!/bin/bash
# LEGACY: sequential all-libraries. Prefer sRNAbench_{LIBRARY}.sbatch (one Slurm job per library).
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name sRNAbench-Sulstoni_newGenome
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
TRIM=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/TrimmedFastq
OUT=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Sulstoni_newGenome
mkdir -p "$OUT"
echo start
java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072570.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR0 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072571.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR1 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072572.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR2 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072573.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR3 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072574.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR4 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072575.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR5 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072576.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR6 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072577.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR7 predict=true species=sulstoniNewGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

echo finished.
