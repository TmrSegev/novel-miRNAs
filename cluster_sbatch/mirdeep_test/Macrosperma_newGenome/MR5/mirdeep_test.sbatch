#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-MR5
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

# Hofstenia-style: run from mirdeep_out/MR5/
miRDeep2.pl /groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma_newGenome/mapper_out/macrosperma_Seq_collapsed_MR5.fasta \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma_newGenome/genome/CMACR.caenorhabditis_macrosperma_JU2083_v2.scaffolds.fna \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma_newGenome/mapper_out/macrosperma_Seq_vs_genome_MR5.arf \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/animalsMature.fa none /groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/animalsHairpin.fa -g -1 \
  2>mirdeep.log
echo finished.
