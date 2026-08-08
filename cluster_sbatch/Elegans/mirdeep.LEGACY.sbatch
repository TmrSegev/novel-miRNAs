#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name mirdeep-Elegans
#SBATCH --output mirdeep-%J.out
#SBATCH --mail-user=stome@post.bgu.ac.il
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
echo "=== CE57 ==="
mkdir -p "$ROOT/mirdeep_out/CE57"
(
  cd "$ROOT/mirdeep_out/CE57"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE57.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE57.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE58 ==="
mkdir -p "$ROOT/mirdeep_out/CE58"
(
  cd "$ROOT/mirdeep_out/CE58"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE58.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE58.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE59 ==="
mkdir -p "$ROOT/mirdeep_out/CE59"
(
  cd "$ROOT/mirdeep_out/CE59"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE59.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE59.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE60 ==="
mkdir -p "$ROOT/mirdeep_out/CE60"
(
  cd "$ROOT/mirdeep_out/CE60"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE60.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE60.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE61 ==="
mkdir -p "$ROOT/mirdeep_out/CE61"
(
  cd "$ROOT/mirdeep_out/CE61"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE61.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE61.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE62 ==="
mkdir -p "$ROOT/mirdeep_out/CE62"
(
  cd "$ROOT/mirdeep_out/CE62"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE62.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE62.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE63 ==="
mkdir -p "$ROOT/mirdeep_out/CE63"
(
  cd "$ROOT/mirdeep_out/CE63"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE63.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE63.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE69 ==="
mkdir -p "$ROOT/mirdeep_out/CE69"
(
  cd "$ROOT/mirdeep_out/CE69"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE69.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE69.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE78 ==="
mkdir -p "$ROOT/mirdeep_out/CE78"
(
  cd "$ROOT/mirdeep_out/CE78"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE78.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE78.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE79 ==="
mkdir -p "$ROOT/mirdeep_out/CE79"
(
  cd "$ROOT/mirdeep_out/CE79"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE79.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE79.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE80 ==="
mkdir -p "$ROOT/mirdeep_out/CE80"
(
  cd "$ROOT/mirdeep_out/CE80"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE80.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE80.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo "=== CE81 ==="
mkdir -p "$ROOT/mirdeep_out/CE81"
(
  cd "$ROOT/mirdeep_out/CE81"
  miRDeep2.pl ../../mapper_out/elegans_Seq_collapsed_CE81.fasta \
    ../../Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
    ../../mapper_out/elegans_Seq_vs_genome_CE81.arf \
    ../../../mirbase_data/animalsMature.fa none ../../../mirbase_data/animalsHairpin.fa -g -1 \
    2>mirdeep.log
)
echo finished.
