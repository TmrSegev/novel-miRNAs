#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-CE63
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

# Hofstenia-style: run from mirdeep_out/CE63/
miRDeep2.pl /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/mapper_out/elegans_Seq_collapsed_CE63.fasta \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/genome/CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/mapper_out/elegans_Seq_vs_genome_CE63.arf \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/animalsMature.fa none /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/animalsHairpin.fa -g -1 \
  2>mirdeep.log
echo finished.
