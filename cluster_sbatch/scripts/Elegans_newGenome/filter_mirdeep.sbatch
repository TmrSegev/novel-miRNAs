#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name filter-mirdeep-Elegans_newGenome
#SBATCH --output filter-mirdeep-%J.out
#SBATCH --mail-user=tzahben@post.bgu.ac.il
#SBATCH --mail-type=ALL
#SBATCH --mem=4G
#SBATCH --cpus-per-task=12

echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

module load anaconda
source activate /home/stome/.conda/envs/my_env

cd "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/mirdeep_out"
for folder in CE57 CE58 CE59 CE60 CE61 CE62 CE63 CE69 CE78 CE79 CE80 CE81; do
  if [ -d "$folder" ]; then
    (cd "$folder" && echo "$folder" && python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepPerLibraryFilter.py -i result*.csv --filter-s 10 --filter-mc 10 --exclude-c 100)
  fi
done
echo finished.
