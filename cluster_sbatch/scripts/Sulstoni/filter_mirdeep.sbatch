#!/bin/bash
#SBATCH --partition main
#SBATCH --time 3-03:30:00
#SBATCH --job-name filter-mirdeep-Sulstoni
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

cd "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/mirdeep_out"
for folder in SR0 SR1 SR2 SR3 SR4 SR5 SR6 SR7; do
  if [ -d "$folder" ]; then
    (cd "$folder" && echo "$folder" && python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepPerLibraryFilter.py -i result*.csv --filter-s 10 --filter-mc 10 --exclude-c 100)
  fi
done
echo finished.
