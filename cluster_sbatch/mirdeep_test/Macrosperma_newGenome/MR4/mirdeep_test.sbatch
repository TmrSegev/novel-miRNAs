#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-MR4
#SBATCH --output mirdeep-%J.out
#SBATCH --mail-user=tzahben@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=8G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

# Hofstenia-style: run from mirdeep_out/MR4/
miRDeep2.pl ../../mapper_out/macrosperma_Seq_collapsed_MR4.fasta \
  ../../genome/CMACR.caenorhabditis_macrosperma_JU2083_v2.scaffolds.fna \
  ../../mapper_out/macrosperma_Seq_vs_genome_MR4.arf \
  ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
  2>mirdeep.log
echo finished.
