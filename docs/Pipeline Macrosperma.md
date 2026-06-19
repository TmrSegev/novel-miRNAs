**Macrosperma**

**Pipeline scripts:** Per-library filter (`srnabenchPerLibraryFilter.py`, `mirdeepPerLibraryFilter.py`; `--filter-mc 10`); unite + GFF in `{base}/Macrosperma/scripts/` (`srnabenchUniteGFF.py`, `mirdeepUniteGFF.py`, `processGoodCandidates.py`). See Pipeline Hofstenia for the two-pass good_candidates workflow. Legacy: `sRNAbenchResultsToGFF3.py`, `mirdeepResultsToGFF3.py`.

Libraries: MR4, MR5, MR6, MR7, MR8 (5 libraries — run miRDeep and sRNAbench **separately per library**, do not combine FASTQs for discovery).

1) **Trimming** **adaptor** in the sequences:

   cutadapt \-a AACTGTAGGCACCATCAAT \--core 2 \-e 0.25 \--discard-untrimmed \-m 17 \-M 26  ../Fastq/SRR13072564.1.fastq \> ../TrimmedFastq/SRR13072564.1\_trimmed.fastq

   

   path:\<basePath\>/Macrosperma/TrimmedFastq/\*.\_trimmed.fastq

   

   Results:

| Library | Reads | Reads with adapters | Reads (after filters) |
| ----- | ----- | ----- | ----- |
| SRR13072564.1 | 69,019,846 | 63,997,300 (92.7%) | 47,186,787 (68.4%) |
| SRR13072565.1 | 65,554,997 | 61,682,578 (94.1%) | 48,478,651 (74.0%) |
| SRR13072566.1 | 57,781,780 | 52,869,020 (91.5%) | 38,837,454 (67.2%) |
| SRR13072567.1 | 61,173,363 | 56,736,384 (92.7%) | 43,294,919 (70.8%) |
| SRR13072568.1 | 76,292,551 | 70,995,798 (93.1%) | 57,286,962 (75.1%) |
| **Total** | **329,822,537** | **306,281,080 (92.8%)** | **235,084,773 (71.2%)** |

   

   

   

2) **genome indexed** by **bowtie-build** \- command:

   bowtie-build \-f CMACR.caenorhabditis\_macrosperma\_JU2083\_v1.scaffolds.fna index/macrospermaGenomeIndexed

   

   path: \<basePath\>/Macrosperma/genome/index/\*

   

3) creating **config.txt**:  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Macrosperma/TrimmedFastq/SRR13072564.1\_trimmed.fastq MR4  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Macrosperma/TrimmedFastq/SRR13072565.1\_trimmed.fastq MR5  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Macrosperma/TrimmedFastq/SRR13072566.1\_trimmed.fastq MR6  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Macrosperma/TrimmedFastq/SRR13072567.1\_trimmed.fastq MR7  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Macrosperma/TrimmedFastq/SRR13072568.1\_trimmed.fastq MR8  
   

path: \<basePath\>/Macrosperma/bash/config.txt

4) **mapper.pl** \- command:

   mapper.pl config.txt \-d \-e \-i \-j \-m \-h \-p ../genome/index/macrospermaGenomeIndexed \-t ../mapper\_out/macrosperma\_Seq\_vs\_genome.arf \-s ../mapper\_out/macrosperma\_Seq\_collapsed.fasta

	  
	parameters:  
	\-e: input file in fastq format.  
	\-d: input file is config.  
	\-p: prefix name of the genome indexed by bowtie.  
	\-i: convert rna to dna (for map againsts genome).  
	\-j: remove sequences thats contain empty value like ‘n’.  
	\-t: print mapping read to the file.  
	\-s: print collapsed reads to the file.  
		\-h: changing the output to fasta format (required).

		Results:

| desc | total | mapped | unmapped | %mapped | %unmapped |
| ----- | ----- | ----- | ----- | ----- | ----- |
| MR4 | 46936385 | 43388384 | 3548001 | 92.441 | 7.559 |
| MR5 | 48241177 | 44455524 | 3785653 | 92.153 | 7.847 |
| MR6 | 38593759 | 35931364 | 2662395 | 93.101 | 6.899 |
| MR7 | 43034731 | 40044631 | 2990100 | 93.052 | 6.948 |
| MR8 | 56925027 | 53050372 | 3874655 | 93.193 | 6.807 |
| **Total** | **233731079** | **216870275** | **16860804** | **92.786** | **7.214** |

		path: \<basePath\>/Macrosperma/mapper\_out/\*

5) **mirDeep2.pl** — run **per library** in `{base}/Macrosperma/mirdeep_out/{MR4|MR5|…}/` (legacy: single combined run below). Filter each folder:

   python mirdeepPerLibraryFilter.py -i result\_\*.csv --filter-s 10 --exclude-c 100 --filter-mc 10

   **Legacy combined run** (superseded):

   

   Results:

   

   

   

| novel miRNAs |  |  | known miRBase miRNAs |  |  |  |  |  |
| :---- | :---- | :---- | :---- | :---- | :---- | ----- | ----- | :---- |
| **miRDeep2 score** | **Predicted by miRDeep2** | **Estimated false positive** | **Estimated true positive** | **In species** | **In data** | **Detected by mirDeep2** | **Estimated signal to noise** | **Excision gearing** |
| 10 | 267 | 44 ± 7 | 223 ± 7 (83 ± 3%) | 37810 | 316 | 289 (91%) | 7 | 1 |
| 9 | 270 | 45 ± 7 | 225 ± 7 (83 ± 3%) | 37810 | 316 | 289 (91%) | 6.9 | 1 |
| 8 | 277 | 46 ± 7 | 231 ± 7 (83 ± 2%) | 37810 | 316 | 290 (92%) | 6.9 | 1 |
| 7 | 286 | 48 ± 7 | 238 ± 7 (83 ± 2%) | 37810 | 316 | 290 (92%) | 6.9 | 1 |
| 6 | 296 | 50 ± 7 | 246 ± 7 (83 ± 2%) | 37810 | 316 | 290 (92%) | 6.9 | 1 |
| 5 | 303 | 52 ± 7 | 251 ± 7 (83 ± 2%) | 37810 | 316 | 290 (92%) | 6.7 | 1 |
| 4 | 321 | 55 ± 8 | 266 ± 8 (83 ± 2%) | 37810 | 316 | 290 (92%) | 6.7 | 1 |
| 3 | 352 | 59 ± 8 | 293 ± 8 (83 ± 2%) | 37810 | 316 | 290 (92%) | 6.7 | 1 |
| 2 | 449 | 73 ± 9 | 376 ± 9 (84 ± 2%) | 37810 | 316 | 298 (94%) | 6.8 | 1 |
| 1 | 567 | 120 ± 10 | 447 ± 10 (79 ± 2%) | 37810 | 316 | 298 (94%) | 5.1 | 1 |
| 0 | 642 | 409 ± 19 | 233 ± 19 (36 ± 3%) | 37810 | 316 | 298 (94%) | 1.7 | 1 |
| \-1 | 732 | 688 ± 24 | 44 ± 23 (6 ± 3%) | 37810 | 316 | 298 (94%) | 1.1 | 1 |
| \-2 | 940 | 943 ± 29 | 10 ± 16 (1 ± 2%) | 37810 | 316 | 298 (94%) | 1 | 1 |
| \-3 | 1652 | 1245 ± 35 | 407 ± 35 (25 ± 2%) | 37810 | 316 | 298 (94%) | 1.4 | 1 |
| \-4 | 2983 | 1934 ± 40 | 1049 ± 40 (35 ± 1%) | 37810 | 316 | 298 (94%) | 1.6 | 1 |
| \-5 | 3729 | 4854 ± 54 | 0 ± 0 (0 ± 0%) | 37810 | 316 | 298 (94%) | 0.8 | 1 |
| \-6 | 4291 | 9257 ± 72 | 0 ± 0 (0 ± 0%) | 37810 | 316 | 298 (94%) | 0.5 | 1 |
| \-7 | 5773 | 12602 ± 83 | 0 ± 0 (0 ± 0%) | 37810 | 316 | 298 (94%) | 0.5 | 1 |
| \-8 | 10053 | 15292 ± 89 | 0 ± 0 (0 ± 0%) | 37810 | 316 | 299 (95%) | 0.7 | 1 |
| \-9 | 13557 | 17253 ± 95 | 0 ± 0 (0 ± 0%) | 37810 | 316 | 299 (95%) | 0.8 | 1 |
| \-10 | 16229 | 18726 ± 95 | 0 ± 0 (0 ± 0%) | 37810 | 316 | 299 (95%) | 0.9 | 1 |

   

   path: \<basePath\>/Macrosperma/mirdeep\_out/\*

   

6) another strategy for prediction \-using sRNAtoolbox,  
   **makeSeqObj.jar** \- preprocess genome \- command:

   java \-jar ../../sRNAtoolboxDB/exec/makeSeqObj.jar ../genome/CMACR.caenorhabditis\_macrosperma\_JU2083\_v1.scaffolds.fna  
     
   path: \<basePath\>/sRNAtoolboxDB/seqOBJ/macrospermaGenomeIndexed.zip  
     
7) **Per-library trimmed FASTQ** — do **not** combine for discovery (legacy `cat SRR* > Macrosperma_final.fastq` superseded).

8) **sRNAbench.jar** — run **per library**. Example MR4:

   java \-jar ../../sRNAtoolboxDB/exec/sRNAbench.jar input=../TrimmedFastq/SRR13072564.1\_trimmed.fastq output=../../sRNAtoolboxDB/out/Macrosperma/Macrosperma\_MR4 predict=true species=macrospermaGenomeIndexed dbPath=/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB hairpin=animalsHairpin.fa mature=animalsMature.fa

   Repeat for MR5–MR8. Filter in each folder:

   python srnabenchPerLibraryFilter.py -i novel.txt -a novel451.txt --filter-mc 10

   **Legacy combined run** (superseded):

   

   Results:

   	Found 458 putative novel

                  Got 445 putative novel after removing redundancies

                  Found 276 putative novel 451-like

                  Got 274 putative 451-like novel after removing redundancies

   

   path: \<basePath\>/Macrosperma/sRNAbench\_out/\*

9) Per-library sRNAbench outputs live under `{base}/sRNAtoolboxDB/out/Macrosperma/Macrosperma_{library}/` (no copy to `Macrosperma/sRNAbench_out/` needed).

**Uniting, good_candidates, and creating GFF3/FASTA**

Working directory: `{base}/Macrosperma/scripts/`

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#  
1. GFF with `--goodcandidates False` → 2. `processGoodCandidates.py` → 3. GFF with `--goodcandidates True` → 4. `compare_genome_to_fasta.py --mode discovery`  
\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

python srnabenchUniteGFF.py -o Macrosperma\_sRNAbench.gff3 -seed ../../mirbase\_data/Seeds.txt --create-fasta Macrosperma\_sRNAbench.fasta -s Macrosperma --goodcandidates False  
python processGoodCandidates.py --tool sRNAbench -s Macrosperma  
python srnabenchUniteGFF.py -o Macrosperma\_sRNAbench.gff3 -seed ../../mirbase\_data/Seeds.txt --create-fasta Macrosperma\_sRNAbench.fasta -s Macrosperma --goodcandidates True  

python mirdeepUniteGFF.py -o Macrosperma\_mirdeep.gff3 --create-fasta Macrosperma\_mirdeep.fasta -seed ../../mirbase\_data/Seeds.txt -s Macrosperma --goodcandidates False  
python processGoodCandidates.py --tool miRDeep -s Macrosperma  
python mirdeepUniteGFF.py -o Macrosperma\_mirdeep.gff3 --create-fasta Macrosperma\_mirdeep.fasta -seed ../../mirbase\_data/Seeds.txt -s Macrosperma --goodcandidates True  

python compare_genome_to_fasta.py --mode discovery --species Macrosperma --dir /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts --genome_fasta /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/genome/CMACR.caenorhabditis\_macrosperma\_JU2083\_v1.scaffolds.fna --gff Macrosperma\_sRNAbench.gff3 --mature Macrosperma\_sRNAbench.fasta --star Macrosperma\_sRNAbench\_star.fasta --hairpin-table sRNAbench\_all\_remaining\_filtered.csv --output sRNAbench\_coord\_check.csv  

paths (final outputs):  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_sRNAbench.gff3  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_mirdeep.gff3  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/sRNAbench\_all\_remaining\_filtered.csv  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/mirdeep\_all\_remaining\_filtered.csv  

	**Filtering (--filter-mc 10):**  
	sRNAbench: max(5pRC,3pRC)\<10 or matureBindings\<14; all novel451 discarded; ncRNA filter; hairpin trim in filter script.  
	miRDeep: `--filter-s 10 --exclude-c 100 --filter-mc 10` per library.  
	Unite step: coordinate overlap dedup + good_candidates (≥2 libraries in 20 bp cluster for Macrosperma).

**For the candidates that are left, we need to mark them as “sense”/”antisense” or “overlap”:**  
“Antisense” miRNAs overlap another miRNA/candidate on the **opposite** strand:   
The one that has higher counts will be marked as “sense” and the other as “antisense” (the precursor will be marked in the pre-miRNA GFF3 files).  
“Overlap” miRNAs overlap another miRNA/candidate on the **same** strand.  
Strands were found by using bedtools intersect, commands:  
Path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/

1. Preprocessing: remove ‘\\t’ from the end of lines in the gffs.  
   sed \-i 's/\\t\*$//' Macrosperma\_mirdeep\_pre\_only.gff3  
   sed \-i 's/\\t\*$//' Macrosperma\_sRNAbench\_pre\_only.gff3  
2. Run intersections.sbatch file. The commands inside:  
   mirdeep-mirdeep bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_mirdeep\_pre\_only.gff3 \-b /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_mirdeep\_pre\_only.gff3 \> miRdeep\_intersect.bed  
     
   sRNAbench-sRNAbench bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_sRNAbench\_pre\_only.gff3 \-b /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_sRNAbench\_pre\_only.gff3 \> sRNAbench\_intersect.bed  
3. Script commands for marking as overlaps or sense/antisense:  
   mirdeep:  
   python overlapSenseAnti.py \--intersections-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/miRdeep\_intersect.bed \--gff /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_mirdeep\_pre\_only.gff3  
     
   sRNAbench:  
   python overlapSenseAnti.py \--intersections-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/sRNAbench\_intersect.bed \--gff /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_sRNAbench\_pre\_only.gff3  
   

The output changes the pre\_only gff files in the respective folders.

**Intersecting miRdeep & sRNAbench & Known Results (bedtools)**  
Finding candidates that are part of the intersection between two tools or known miRNAs.  
path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/  
All commands documented in \<path\>/Command.txt

1. When applying intersect force use “-s” for strandness, \-f for minimum overlap.  
   f paramater was chosen based on the “cleanliness” of the canditate. In miRGeneDB and miRDeep the candidates are pures. In sRNAbench, each candidate has an extra head and tail of 11 bases each. miRBase has extra head and tail, sometimes longer.  
   f=0.6 was chosen to allow for the difference between sRNAbench and a pure candidate of roughly 50 bases (50/72=0.694). f=0.4 was chosen to allow for the differences between miRBase and the other databases.  
     
   Apply intersect on:  
   

| sRNAbench | miRDeep | f 0.6 |
| :---- | :---- | :---- |

   

| miRDeep | sRNAbench | f 0.6 |
| :---- | :---- | :---- |

**STAR Alignment**

1)   **Genome indexing** by **STAR** \- command:

   STAR \--runMode genomeGenerate \--runThreadN 16 \--genomeDir ../STAR/genome\_index/ \--genomeSAindexNbases 11  \--genomeFastaFiles ../Genome/CMACR.caenorhabditis\_macrosperma\_JU2083\_v1.scaffolds.fna

   

   path: \<basePath\>/Macrosperma/STAR/genome\_index/\*

   

2)  **Align** the combined file (including all Macrosperma libraries) **to genome** by **STAR** \- command for the first library:

   STAR \--genomeDir ../STAR/genome\_index/ \--readFilesIn ../TrimmedFastq/SRR13072557.1\_trimmed.fastq \--outFileNamePrefix ../STAR/align\_to\_genome/CE57/Macrosperma\_ \--outFilterMultimapNmax 20 \--runThreadN 16 \--outSAMtype SAM

   

3) **featureCounts** for **sRNAbench gff3 & miRDeep gff3** \- command:  
   Run the “featurecounts\_\<name\>\_sep.sbatch” file in path.  
   sRNAbench:

   featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_sRNAbench.gff3 \-o ../counts\_sep/miRNA\_sRNAbench\_counts.txt ../STAR/align\_to\_genome/MR8/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR7/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR6/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR5/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR4/Macrosperma\_Aligned.out.sam  
   

	mirdeep:

featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_mirdeep.gff3 \-o ../counts\_sep/miRNA\_miRdeep\_counts.txt ../STAR/align\_to\_genome/MR8/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR7/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR6/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR5/Macrosperma\_Aligned.out.sam ../STAR/align\_to\_genome/MR4/Macrosperma\_Aligned.out.sam

**Flanked precursor counts (m/pre ratio):**

python add_flank_to_GFF.py -s Macrosperma

featureCounts \-F GFF \-t pre\_miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_mirdeep\_flanked\_pre.gff3 \-o ../counts\_sep/miRNA\_miRdeep\_counts\_flanked.txt \<MR4–MR8 SAM files\>

featureCounts \-F GFF \-t pre\_miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_sRNAbench\_flanked\_pre.gff3 \-o ../counts\_sep/miRNA\_sRNAbench\_counts\_flanked.txt \<MR4–MR8 SAM files\>

**BLAST** 

1. Creating a DB of known miRNAs of nematodes only.  
   1. Create blast DB, command:  
      makeblastdb \-in ../BLAST\_DB/Caenorhabditis\_pre\_miRNA.fasta \-title miRNADB \-dbtype nucl \-out ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB  
2. Blast mature results from miRdeep and sRNAbench. Commands:  
   blastn \-query ../../Charles\_seq/Macrosperma/scripts/Macrosperma\_mirdeep.fasta \-db ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB \-out ../queries/Macrosperma/miRdeep\_blastn\_compact \-outfmt 6 \-evalue 10 \-task blastn-short  
   blastn \-query /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/Macrosperma\_sRNAbench.fasta \-db ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB \-out ../queries/Macrosperma/sRNAbench\_blastn\_compact \-outfmt 6 \-evalue 10 \-task blastn-short  
   path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/bash/  
   

**Generate Intersections Table**

1. Path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/  
2. Command for script:  
   python ‏‏‏‏intersectionsTable.py \-s macrosperma \--mirdeep-inter-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/miRdeep\_sRNAbench\_intersect.bed \--sRNAbench-inter-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/sRNAbench\_miRdeep\_intersect.bed \--blast-mirdeep /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/queries/Macrosperma/miRdeep\_blastn\_compact \--blast-sRNAbench /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/queries/Macrosperma/sRNAbench\_blastn\_compact \--fc-mirdeep /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/counts\_sep/miRNA\_miRdeep\_counts.txt \--fc-pre-mirdeep /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/counts\_sep/miRNA\_miRdeep\_counts\_flanked.txt \--fc-sRNAbench /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/counts\_sep/miRNA\_sRNAbench\_counts.txt \--fc-pre-sRNAbench /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/counts\_sep/miRNA\_sRNAbench\_counts\_flanked.txt \-rm /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/mirdeep\_all\_remaining\_filtered.csv \-rs /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/scripts/sRNAbench\_all\_remaining\_filtered.csv \-l MR8,MR7,MR6,MR5,MR4 \--sum-fc-thres 100  
   Output:  
   intersections\_table\_script.xlsx  
   Merges miRDeep results with blast, featurecounts, sRNAbench  
   Merges sRNAbench results with blast, featurecounts, miRdeep

**Generate all candidates fasta**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/

command:

python allCandidatesFasta.py \--all /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/intersections\_table\_macrosperma.xlsx \-s Macrosperma

**Feature Engineering with Ziv’s code**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Ziv\_Features

command:

python Ziv\_feature\_SOS.py \--precursors /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/all\_candidates\_hairpin.fasta \--mature /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/all\_candidates\_mature.fasta \--star /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/all\_candidates\_star.fasta \--species Macrosperma \--all-remaining /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/intersections\_table\_macrosperma.xlsx

**Statistical Analysis**

Path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Macrosperma/

command:

python statistics.py \--all /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Macrosperma.xlsx \-s Macrosperma

