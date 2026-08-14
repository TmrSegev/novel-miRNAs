#!/bin/bash

#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name cutadapt-macrosperma
#SBATCH --output cutadapt-%J.out
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

# Same reads as the old-genome track; mapper/STAR/sRNAbench reuse this TrimmedFastq.
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/Fastq
OUT=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq
mkdir -p "$OUT"

echo start
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072564.1.fastq > ${OUT}/SRR13072564.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072565.1.fastq > ${OUT}/SRR13072565.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072566.1.fastq > ${OUT}/SRR13072566.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072567.1.fastq > ${OUT}/SRR13072567.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072568.1.fastq > ${OUT}/SRR13072568.1_trimmed.fastq
echo finished.
