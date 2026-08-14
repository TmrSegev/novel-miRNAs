#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-Sulstoni_newGenome
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
echo "=== SR0 ==="
mkdir -p "$ROOT/mirdeep_out/SR0"
(
  cd "$ROOT/mirdeep_out/SR0"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR0.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR0.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== SR1 ==="
mkdir -p "$ROOT/mirdeep_out/SR1"
(
  cd "$ROOT/mirdeep_out/SR1"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR1.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR1.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== SR2 ==="
mkdir -p "$ROOT/mirdeep_out/SR2"
(
  cd "$ROOT/mirdeep_out/SR2"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR2.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR2.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== SR3 ==="
mkdir -p "$ROOT/mirdeep_out/SR3"
(
  cd "$ROOT/mirdeep_out/SR3"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR3.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR3.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== SR4 ==="
mkdir -p "$ROOT/mirdeep_out/SR4"
(
  cd "$ROOT/mirdeep_out/SR4"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR4.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR4.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== SR5 ==="
mkdir -p "$ROOT/mirdeep_out/SR5"
(
  cd "$ROOT/mirdeep_out/SR5"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR5.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR5.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== SR6 ==="
mkdir -p "$ROOT/mirdeep_out/SR6"
(
  cd "$ROOT/mirdeep_out/SR6"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR6.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR6.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== SR7 ==="
mkdir -p "$ROOT/mirdeep_out/SR7"
(
  cd "$ROOT/mirdeep_out/SR7"
  miRDeep2.pl ../../mapper_out/sulstoni_Seq_collapsed_SR7.fasta \
    ../../genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna \
    ../../mapper_out/sulstoni_Seq_vs_genome_SR7.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo finished.
