#!/bin/bash
### sbatch config parameters must start with #SBATCH and must precede any other command. to ignore just add another # - like ##SBATCH
#SBATCH --partition main ### partition name where to run a job. Use ‘main’ unless qos is required. qos partitions ‘rtx3090’ ‘rtx2080’ ‘gtx1080’
#SBATCH --time 0-19:00:00 ### limit the time of job running. Make sure it is not greater than the partition time limit (7 days)!! Format: D-H:MM:SS
#SBATCH --job-name mirge ### name of the job. replace my_job with your desired job name
#SBATCH --output mirge-id-%J.out ### output log for running job - %J is the job number variable
##SBATCH --mail-user=stome@post.bgu.ac.il ### user’s email for sending job status notifications
#SBATCH --mail-type=BEGIN,END,FAIL ### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --gpus=0 ### number of GPUs. Choosing type e.g.: #SBATCH --gpus=gtx_1080:1 , or rtx_2080, or rtx_3090 . Allocating more than 1 requires the IT team’s permission
##SBATCH --tasks=1 # 1 process – use for processing of few programs concurrently in a job (with srun). Use just 1 otherwise

echo "Job start"
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/mirge_env				### activate a conda environment, replace my_env with your conda environmen
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/
miRge3.0 \
  -s /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/GA3.filtered.fastq \
  -lib /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/ \
  -on Hofstenia \
  -db miRBase \
  -o /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/results \
  -gff \
  -bam \
  -ie \
  -cpu 8

echo "Job finished"