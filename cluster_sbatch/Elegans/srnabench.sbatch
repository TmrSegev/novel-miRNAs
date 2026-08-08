#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name sRNAbench-Elegans
#SBATCH --output sRNAbench-%J.out
#SBATCH --mail-user=stome@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda

source activate /home/stome/.conda/envs/my_env

DB=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB
TRIM=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/TrimmedFastq
OUT=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Elegans
mkdir -p "$OUT"
echo start
java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072557.1_trimmed.fastq \
  output=${OUT}/Elegans_CE57 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072558.1_trimmed.fastq \
  output=${OUT}/Elegans_CE58 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072559.1_trimmed.fastq \
  output=${OUT}/Elegans_CE59 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072560.1_trimmed.fastq \
  output=${OUT}/Elegans_CE60 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072561.1_trimmed.fastq \
  output=${OUT}/Elegans_CE61 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072562.1_trimmed.fastq \
  output=${OUT}/Elegans_CE62 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072563.1_trimmed.fastq \
  output=${OUT}/Elegans_CE63 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072569.1_trimmed.fastq \
  output=${OUT}/Elegans_CE69 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072578.1_trimmed.fastq \
  output=${OUT}/Elegans_CE78 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072579.1_trimmed.fastq \
  output=${OUT}/Elegans_CE79 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072580.1_trimmed.fastq \
  output=${OUT}/Elegans_CE80 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=${TRIM}/SRR13072581.1_trimmed.fastq \
  output=${OUT}/Elegans_CE81 predict=true species=elegansGenomeIndexed \
  dbPath=${DB} hairpin=animalsHairpin.fa mature=animalsMature.fa

echo finished.
