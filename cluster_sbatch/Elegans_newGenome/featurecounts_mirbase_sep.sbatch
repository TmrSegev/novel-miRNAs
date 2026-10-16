#!/bin/bash

###sbatch configuration parameters must start with #SBATCH and must precede any other commands.
### To ignore, just add another # - like so: ##SBATCH
################################################################################################

#SBATCH --partition main			### specify partition name where to run a job. short: 7 days limit; gtx1080: 7 days; debug: 2 hours limit and 1 job at a time
#SBATCH --time 3-03:30:00			### limit the time of job running. Make sure it is not greater than the partition time limit!! Format: D-H:MM:SS
#SBATCH --job-name featurecounts_mirbase			### name of the job
#SBATCH --output featurecounts_sep_mirbase-%J.out			### output log for running job - %J for job number
##SBATCH --mail-user=stome@post.bgu.ac.il	### user's email for sending job status messages
##SBATCH --mail-type=ALL			### conditions for sending the email. ALL,BEGIN,END,FAIL, REQUEU, NONE
#SBATCH --mem=32G				### ammount of RAM memory
#SBATCH --cpus-per-task=12			### number of CPU cores

### Print some data to output file ###
echo `date`
echo -e "\nSLURM_JOBID:\t\t" $SLURM_JOBID
echo -e "SLURM_JOB_NODELIST:\t" $SLURM_JOB_NODELIST "\n\n"

### Start your code below ####
module load anaconda				### load anaconda module (must be present when working with conda environments)
source activate /home/stome/.conda/envs/my_env				### activate a conda environment, replace my_env with your conda environmen

mkdir -p "/groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/counts_sep"
echo start ## execute python script – replace with your own command 
set -e
# Same featureCounts call as miRDeep/sRNAbench/Hofstenia: -t miRNA -g ID, no -F.
# Those GFF3 miRNA rows are only "ID=<id>". cel_mirbase_seq.gff3 appends a raw
# sequence and 5p/3p with no key, and featureCounts 2.0.1 then cannot read ID=.
# Rewrite miRNA rows to that GFF3 shape. Phase 10 -mgff stays the original file.
GFF_IN=/groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/cel_mirbase_seq.gff3
GFF_OUT=/groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/counts_sep/cel_mirbase_miRNA.gff3
awk -F'\t' 'BEGIN{OFS="\t"; print "##gff-version 3"}
$0 !~ /^#/ && $3 == "miRNA" && NF >= 9 {
  id = $9
  sub(/;.*/, "", id)
  if (id !~ /^ID=/) id = "ID=" id
  if (id != "ID=") print $1, $2, "miRNA", $4, $5, $6, $7, $8, id
}' "$GFF_IN" > "$GFF_OUT"
nfeat=$(awk 'END{print NR-1}' "$GFF_OUT")
echo "GFF3 miRNA features: $nfeat"
[[ "$nfeat" -gt 0 ]]
featureCounts -t miRNA -g ID -O -s 1 -M \
  -a "$GFF_OUT" \
  -o /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/counts_sep/miRNA_mirbase_counts.txt \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE81/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE80/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE69/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE63/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE62/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE61/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE60/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE59/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE58/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE57/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE79/Elegans_Aligned.out.sam \
  /groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans_newGenome/STAR/align_to_genome/CE78/Elegans_Aligned.out.sam
echo finished.
