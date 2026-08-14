#!/bin/bash

#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name cutadapt-sulstoni
#SBATCH --output cutadapt-%J.out
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

# Same reads as the old-genome track; mapper/STAR/sRNAbench reuse this TrimmedFastq.
FQ=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/Fastq
OUT=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/TrimmedFastq
mkdir -p "$OUT"

echo start
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072570.1.fastq > ${OUT}/SRR13072570.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072571.1.fastq > ${OUT}/SRR13072571.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072572.1.fastq > ${OUT}/SRR13072572.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072573.1.fastq > ${OUT}/SRR13072573.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072574.1.fastq > ${OUT}/SRR13072574.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072575.1.fastq > ${OUT}/SRR13072575.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072576.1.fastq > ${OUT}/SRR13072576.1_trimmed.fastq
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 ${FQ}/SRR13072577.1.fastq > ${OUT}/SRR13072577.1_trimmed.fastq
echo finished.
