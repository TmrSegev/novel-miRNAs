#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name filter-sRNAbench-Sulstoni_newGenome
#SBATCH --output filter-sRNAbench-%J.out
#SBATCH --mail-user=tzahben@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

cd "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Sulstoni_newGenome"
for folder in Sulstoni_SR0 Sulstoni_SR1 Sulstoni_SR2 Sulstoni_SR3 Sulstoni_SR4 Sulstoni_SR5 Sulstoni_SR6 Sulstoni_SR7; do
  if [ -d "$folder" ]; then
    (cd "$folder" && echo "$folder" && python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py -i novel.txt -a novel451.txt --filter-mc 10)
  fi
done
echo finished.
