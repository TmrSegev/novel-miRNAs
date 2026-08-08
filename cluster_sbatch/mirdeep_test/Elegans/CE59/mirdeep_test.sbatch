#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-CE59
#SBATCH --output mirdeep-%J.out
#SBATCH --mail-user=stome@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=8G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

# Hofstenia-style: run from mirdeep_out/CE59/
miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE59.fasta \
  ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
  ../../mapper_out/elegans_Seq_vs_genome_CE59.arf \
  ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
  2>mirdeep.log
echo finished.
