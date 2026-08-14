#!/bin/bash

#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name cutadapt-elegans
#SBATCH --output cutadapt-%J.out
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

# Same reads as the old-genome track; mapper/STAR/sRNAbench reuse this TrimmedFastq.
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/Fastq
OUT=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/TrimmedFastq
mkdir -p "$OUT"

echo start
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072557.1.fastq > ${OUT}/SRR13072557.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072558.1.fastq > ${OUT}/SRR13072558.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072559.1.fastq > ${OUT}/SRR13072559.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072560.1.fastq > ${OUT}/SRR13072560.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072561.1.fastq > ${OUT}/SRR13072561.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072562.1.fastq > ${OUT}/SRR13072562.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072563.1.fastq > ${OUT}/SRR13072563.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072569.1.fastq > ${OUT}/SRR13072569.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072578.1.fastq > ${OUT}/SRR13072578.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072579.1.fastq > ${OUT}/SRR13072579.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072580.1.fastq > ${OUT}/SRR13072580.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 12 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072581.1.fastq > ${OUT}/SRR13072581.1_trimmed.fastq
echo finished.
