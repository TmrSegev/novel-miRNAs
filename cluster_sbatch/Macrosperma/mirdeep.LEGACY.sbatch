#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-Macrosperma
#SBATCH --output mirdeep-%J.out
#SBATCH --mail-user=tzahben@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=128G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda

source activate /home/stome/.conda/envs/my_env

# Per-library miRDeep2 (Hofstenia-style). Requires mapper_out/*_{LIBRARY}.*
ROOT="$(cd .. && pwd)"
echo start
echo "=== MR4 ==="
mkdir -p "$ROOT/mirdeep_out/MR4"
(
  cd "$ROOT/mirdeep_out/MR4"
  miRDeep2.pl ../../mapper_out/macrosperma_Seq_collapsed_MR4.fasta \
    ../../genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna \
    ../../mapper_out/macrosperma_Seq_vs_genome_MR4.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== MR5 ==="
mkdir -p "$ROOT/mirdeep_out/MR5"
(
  cd "$ROOT/mirdeep_out/MR5"
  miRDeep2.pl ../../mapper_out/macrosperma_Seq_collapsed_MR5.fasta \
    ../../genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna \
    ../../mapper_out/macrosperma_Seq_vs_genome_MR5.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== MR6 ==="
mkdir -p "$ROOT/mirdeep_out/MR6"
(
  cd "$ROOT/mirdeep_out/MR6"
  miRDeep2.pl ../../mapper_out/macrosperma_Seq_collapsed_MR6.fasta \
    ../../genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna \
    ../../mapper_out/macrosperma_Seq_vs_genome_MR6.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== MR7 ==="
mkdir -p "$ROOT/mirdeep_out/MR7"
(
  cd "$ROOT/mirdeep_out/MR7"
  miRDeep2.pl ../../mapper_out/macrosperma_Seq_collapsed_MR7.fasta \
    ../../genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna \
    ../../mapper_out/macrosperma_Seq_vs_genome_MR7.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== MR8 ==="
mkdir -p "$ROOT/mirdeep_out/MR8"
(
  cd "$ROOT/mirdeep_out/MR8"
  miRDeep2.pl ../../mapper_out/macrosperma_Seq_collapsed_MR8.fasta \
    ../../genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna \
    ../../mapper_out/macrosperma_Seq_vs_genome_MR8.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo finished.
