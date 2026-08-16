#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name sRNAbench-SR5
#SBATCH --output sRNAbench-SR5-%J.out
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
OUT=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Sulstoni
mkdir -p "$OUT"
echo start
java -jar /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072575.1_trimmed.fastq \
  output=${OUT}/Sulstoni_SR5 predict=true species=sulstoniGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa
echo finished.
