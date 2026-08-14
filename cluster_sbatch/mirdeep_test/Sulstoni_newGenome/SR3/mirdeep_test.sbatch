#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-SR3
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

# Hofstenia-style: run from mirdeep_out/SR3/
miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR3.fasta \
  ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
  ../../mapper_out/sulstoni_Seq_vs_genome_SR3.arf \
  ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
  2>mirdeep.log
echo finished.
