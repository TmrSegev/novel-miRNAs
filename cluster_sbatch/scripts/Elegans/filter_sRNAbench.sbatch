#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name filter-sRNAbench-Elegans
#SBATCH --output filter-sRNAbench-%J.out
#SBATCH --mail-user=stome@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

cd "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Elegans"
for folder in Elegans_CE57 Elegans_CE58 Elegans_CE59 Elegans_CE60 Elegans_CE61 Elegans_CE62 Elegans_CE63 Elegans_CE69 Elegans_CE78 Elegans_CE79 Elegans_CE80 Elegans_CE81; do
  if [ -d "$folder" ]; then
    (cd "$folder" && echo "$folder" && python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py -i novel.txt -a novel451.txt --filter-mc 10)
  fi
done
echo finished.
