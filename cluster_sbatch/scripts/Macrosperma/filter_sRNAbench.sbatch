#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name filter-sRNAbench-Macrosperma
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

cd "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Macrosperma"
for folder in Macrosperma_MR4 Macrosperma_MR5 Macrosperma_MR6 Macrosperma_MR7 Macrosperma_MR8; do
  if [ -d "$folder" ]; then
    (cd "$folder" && echo "$folder" && python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py -i novel.txt -a novel451.txt --filter-mc 10)
  fi
done
echo finished.
