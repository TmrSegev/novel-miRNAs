**Sulstoni**

Check position in genome:  
samtools faidx ../genome/CSULS.caenorhabditis\_sulstoni\_JU2788\_v1.scaffolds.fna CSULS.scaffold02010:60847-60868

basePath \= /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq

creating gff3 from sRNAbench output: [https://github.com/tzahy4530/sRNAbenchTableToGFF3](https://github.com/tzahy4530/sRNAbenchTableToGFF3)  
creating gff3 from miRDeep2 output: [https://github.com/tzahy4530/miRDeepResultsToGFF3](https://github.com/tzahy4530/miRDeepResultsToGFF3)  
classify fastqc overpresented sequences: [https://github.com/tzahy4530/fastqcDataClassify](https://github.com/tzahy4530/fastqcDataClassify)

filter fasta files with only matures \= sed \-e '/|s|/,+1d' Sulstoni\_mirdeep.txt \> Sulstoni\_mirdeep\_mature.fasta

Sulstoni:

1) **Trimming** **adaptor** in the sequences:

   cutadapt \-a AACTGTAGGCACCATCAAT \--core 2 \-e 0.25 \--discard-untrimmed \-m 17 \-M 26 ../Fastq/SRR13072570.[https://github.com/tzahy4530/fastqcDataClassify](https://github.com/tzahy4530/fastqcDataClassify)1.fastq \> ../TrimmedFastq/SRR13072570.1\_trimmed.fastq

   

   path:\<basePath\>/Sulstoni/TrimmedFastq/\*.\_trimmed.fastq

   

   Results:

| Library | Reads | Reads with adapters | Reads (after filters) |
| ----- | ----- | ----- | ----- |
| SRR13072570.1 | 18,255,026 | 16,483,389 (90.3%) | 12,910,888 (70.7%) |
| SRR13072571.1 | 16,365,766 | 14,287,310 (87.3%) | 10,645,401 (65.0%) |
| SRR13072572.1 | 15,745,910 | 13,779,371 (87.5%) | 10,349,698 (65.7%) |
| SRR13072573.1 | 16,661,940 | 15,018,391 (90.1%) | 11,333,566 (68.0%) |
| SRR13072574.1 | 15,555,055 | 13,855,319 (89.1%) | 10,482,433 (67.4%) |
| SRR13072575.1 | 17,673,206 | 16,313,109 (92.3%) | 12,723,653 (72.0%) |
| SRR13072576.1 | 21,300,249 | 19,858,596 (93.2%) | 15,388,151 (72.2%) |
| SRR13072577.1 | 20,601,843 | 18,775,839 (91.1%) | 13,654,108 (66.3%) |
| **Total** | **142,158,995** | **128,371,324 (90.3%)** | **97,487,898 (68.5%)** |

   

   

   

2) **genome indexed** by **bowtie-build** \- command:

   bowtie-build \-f CSULS.caenorhabditis\_sulstoni\_JU2788\_v1.scaffolds.fna index/sulstoniGenomeIndexed

   

   path: \<basePath\>/Sulstoni/genome/index/\*

   

3) creating **config.txt**:  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072570.1\_trimmed.fastq SR0  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072571.1\_trimmed.fastq SR1  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072572.1\_trimmed.fastq SR2  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072573.1\_trimmed.fastq SR3  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072574.1\_trimmed.fastq SR4  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072575.1\_trimmed.fastq SR5  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072576.1\_trimmed.fastq SR6  
   /storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/TrimmedFastq/SRR13072577.1\_trimmed.fastq SR7  
     
4) **mapper.pl** \- command:

   mapper.pl config.txt \-d \-e \-i \-j \-m \-h \-p ../genome/index/sulstoniGenomeIndexed \-t ../mapper\_out/Sulstoni\_Seq\_vs\_genome.arf \-s ../mapper\_out/Sulstoni\_Seq\_collapsed.fasta

	  
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
| SR0 | 12812678 | 11879397 | 933281 | 92.716 | 7.284 |
| SR1 | 10558379 | 9800573 | 757806 | 92.823 | 7.177 |
| SR2 | 10266516 | 9555650 | 710866 | 93.076 | 6.924 |
| SR3 | 11241612 | 10477088 | 764524 | 93.199 | 6.801 |
| SR4 | 10393736 | 9707372 | 686364 | 93.396 | 6.604 |
| SR5 | 12596537 | 11734049 | 862488 | 93.153 | 6.847 |
| SR6 | 15238335 | 14199038 | 1039297 | 93.180 | 6.820 |
| SR7 | 13501053 | 12519910 | 981143 | 92.733 | 7.267 |
| **Total** | **96608846** | **89873077** | **6735769** | **93.028** | **6.972** |

		path: \<basePath\>/Sulstoni/mapper\_out/\*

5) **mirDeep2.pl** \- allSeq\_config with Hairpain & Mature \- command:

   miRDeep2.pl ../mapper\_out/Sulstoni\_Seq\_collapsed.fasta ../genome/CSULS.caenorhabditis\_sulstoni\_JU2788\_v1.scaffolds.fna ../mapper\_out/Sulstoni\_Seq\_vs\_genome.arf ../../mirbase\_data/animalsMature.fa none ../../mirbase\_data/animalsHairpin.fa \-g \-1 2\>10.7.2021.mirdeep.log

   

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

   

path: \<basePath\>/Sulstoni/mirdeep\_out/\*

6) another strategy for prediction \-using sRNAtoolbox,  
   **makeSeqObj.jar** \- preprocess genome \- command:

   java \-jar ../../sRNAtoolboxDB/exec/makeSeqObj.jar ../genome/CSULS.caenorhabditis\_sulstoni\_JU2788\_v1.scaffolds.fna  
     
   path: \<basePath\>/sRNAtoolboxDB/seqOBJ/sulstoniGenomeIndexed.zip  
     
7) **Combining** the trimmed **Sulstoni libraries** \- command:  
   	cat SRR\* \> Sulstoni\_final.fastq  
   	

   path: \<basePath\>/Sulstoni/TrimmedFastq/Sulstoni\_final.fastq

     
8) **sRNAbench.jar** \- align to genome and prediction for each library \- command:

   java \-jar ../../sRNAtoolboxDB/exec/sRNAbench.jar input=../TrimmedFastq/Sulstoni\_final.fastq  output=../../sRNAtoolboxDB/out/Sulstoni predict=true species=sulstoniGenomeIndexed dbPath=/storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB hairpin=animalsHairpin.fa mature=animalsMature.fa

   

   Results:

   	Found 301 putative novel

                  Got 292 putative novel after removing redundancies

                  Found 39 putative novel 451-like

                  Got 37 putative 451-like novel after removing redundancies

   

   path: \<basePath\>/Sulstoni/sRNAbench\_out/\*

9) Copy folder from sRNAtoolboxDB/out/Sulstoni/ to Charles\_seq/Sulstoni/sRNAbench\_out/:  
   Commands:  
   cd /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB/out/  
     
   cp \-r ./Sulstoni/\* /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/sRNAbench\_out/  
   

**Creating gff3 & fasta from sRNAbench & mirdeep Output**

1) sRNAbench command:

   python sRNAbenchResultsToGFF3.py \-i novel.txt \-a novel451.txt \-o Sulstoni\_sRNAbench.gff3 \-seed ../../mirbase\_data/Seeds.txt  \--create-fasta Sulstoni\_sRNAbench.fasta \--filter-mc 100 \-s Sulstoni

   miRDeep command:

   python mirdeepResultsToGFF3.py \-i result\_29\_10\_2021\_t\_12\_09\_21.csv \--filter-s 10 \--exclude-c 1000 \--filter-mc 100 \-o Sulstoni\_mirdeep.gff3 \--create-fasta Sulstoni\_mirdeep.fasta \-seed ../../mirbase\_data/Seeds.txt \-s Sulstoni

	  
paths:  
\<basePath\>/Sulstoni/sRNAbench\_out/  
\<basePath\>/Sulstoni/mirdeep\_out/

**For the candidates that are left, we need to mark them as “sense”/”antisense” or “overlap”:**  
“Antisense” miRNAs overlap another miRNA/candidate on the opposite strand:   
The one that has higher counts will be marked as “sense” and the other as “antisense” (the precursor will be marked in the pre-miRNA GFF3 files).  
“Overlap” miRNAs overlap another miRNA/candidate on the **same** strand.  
“Sense” and “antisense” strands were found by using bedtools intersect, commands:

1. Preprocessing: remove ‘\\t’ from the end of lines in the gffs.  
   sed \-i 's/\\t\*$//' Sulstoni\_mirdeep\_pre\_only.gff3  
   sed \-i 's/\\t\*$//' Sulstoni\_sRNAbench\_pre\_only.gff3  
1. Run intersections.sbatch file. The commands inside:  
   mirdeep-mirdeep bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/mirdeep\_out/Sulstoni\_mirdeep\_pre\_only.gff3 \-b /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/mirdeep\_out/Sulstoni\_mirdeep\_pre\_only.gff3 \> miRdeep\_intersect.bed  
     
   sRNAbench-sRNAbench bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/sRNAbench\_out/Sulstoni\_sRNAbench\_pre\_only.gff3 \-b /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/sRNAbench\_out/Sulstoni\_sRNAbench\_pre\_only.gff3 \> sRNAbench\_intersect.bed  
2. Script commands for marking as overlaps or sense/antisense:  
   mirdeep:  
   python overlapSenseAnti.py \--intersections-table /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/miRdeep\_intersect.bed \--gff /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/mirdeep\_out/Sulstoni\_mirdeep\_pre\_only.gff3  
     
   sRNAbench:  
   python overlapSenseAnti.py \--intersections-table /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/sRNAbench\_intersect.bed \--gff /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/sRNAbench\_out/Sulstoni\_sRNAbench\_pre\_only.gff3  
   

The output changes the pre\_only gff files in the respective folders.

**Intersecting miRdeep & sRNAbench & Known Results (bedtools)**  
Finding candidates that are part of the intersection between two tools or known miRNAs.  
path: /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/  
All commands documented in \<path\>/Command.txt

2. When applying intersect force use “-s” for strandness, \-f for minimum overlap.  
   f paramater was chosen based on the “cleanliness” of the canditate. In miRGeneDB and miRDeep the candidates are pures. In sRNAbench, each candidate has an extra head and tail of 11 bases each. miRBase has extra head and tail, sometimes longer.  
   f=0.6 was chosen to allow for the difference between sRNAbench and a pure candidate of roughly 50 bases (50/72=0.694). f=0.5 was chosen to allow for the differences between miRBase and the other databases.  
     
   Apply intersect on:  
   

| sRNAbench | miRDeep | f 0.6 |
| :---- | :---- | :---- |

   

| miRDeep | sRNAbench | f 0.6 |
| :---- | :---- | :---- |

**STAR Alignment**

1)   **Genome indexing** by **STAR** \- command:

   STAR \--runMode genomeGenerate \--runThreadN 16 \--genomeDir ../STAR/genome\_index/ \--genomeSAindexNbases 11  \--genomeFastaFiles ../Genome/CSULS.caenorhabditis\_sulstoni\_JU2788\_v1.scaffolds.fna

   

   path: \<basePath\>/Elegans/STAR/genome\_index/\*

   

2)  **Align** the combined file (including all Elegans libraries) **to genome** by **STAR** \- command for the first library:

   STAR \--genomeDir ../STAR/genome\_index/ \--readFilesIn ../TrimmedFastq/SRR13072557.1\_trimmed.fastq \--outFileNamePrefix ../STAR/align\_to\_genome/CE57/Sulstoni\_ \--outFilterMultimapNmax 20 \--runThreadN 16 \--outSAMtype SAM

   

3) **featureCounts** for **sRNAbench gff3 & miRDeep gff3** \- command:  
   mirdeep:

   featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/mirdeep\_out/Sulstoni\_mirdeep.gff3 \-o ../counts\_sep/miRNA\_miRdeep\_counts.txt ../STAR/align\_to\_genome/SR7/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR6/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR5/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR4/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR3/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR2/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR1/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR0/Sulstoni\_Aligned.out.sam  
   

	sRNAbench:  
		featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/sRNAbench\_out/Sulstoni\_sRNAbench.gff3 \-o ../counts\_sep/miRNA\_sRNAbench\_counts.txt ../STAR/align\_to\_genome/SR7/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR6/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR5/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR4/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR3/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR2/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR1/Sulstoni\_Aligned.out.sam ../STAR/align\_to\_genome/SR0/Sulstoni\_Aligned.out.sam

**BLAST** 

1. Creating a DB of known miRNAs of nematodes only.  
   1. Create blast DB, command:  
      makeblastdb \-in ../BLAST\_DB/Caenorhabditis\_pre\_miRNA.fasta \-title miRNADB \-dbtype nucl \-out ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB  
2. Blast mature results from miRdeep and sRNAbench. Commands:  
   blastn \-query /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/mirdeep\_out/Sulstoni\_mirdeep.fasta \-db ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB \-out ../queries/Sulstoni/miRdeep\_blastn\_compact \-outfmt 6 \-evalue 10 \-task blastn-short  
   blastn \-query /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/sRNAbench\_out/Sulstoni\_sRNAbench.fasta \-db ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB \-out ../queries/Sulstoni/sRNAbench\_blastn\_compact \-outfmt 6 \-evalue 10 \-task blastn-short  
   path: /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/bash/blast\_sulstoni\_queries.sbatch  
   

**Generate Intersections Table**

1. Path: /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni  
2. Command for script:  
   python ‏‏‏‏intersectionsTable.py \-s sulstoni \--mirdeep-inter-table /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/miRdeep\_sRNAbench\_intersect.bed \--sRNAbench-inter-table /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/sRNAbench\_miRdeep\_intersect.bed \--blast-mirdeep /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/queries/Sulstoni/miRdeep\_blastn\_compact \--blast-sRNAbench /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/queries/Sulstoni/sRNAbench\_blastn\_compact \--fc-mirdeep /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/counts\_sep/miRNA\_miRdeep\_counts.txt \--fc-sRNAbench /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/counts\_sep/miRNA\_sRNAbench\_counts.txt \-r1m /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/mirdeep\_out/remaining\_file\_1.csv \-r2m /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/mirdeep\_out/remaining\_file\_2.csv \-rs /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Sulstoni/sRNAbench\_out/sRNAbench\_remaining.csv \-l SR7,SR6,SR5,SR4,SR3,SR2,SR1,SR0 –sum-fc-thres 100  
   Output:  
   intersections\_table\_sulstoni.xlsx  
   Merges miRDeep results with blast, featurecounts, sRNAbench, miRbase and miRGeneDB.  
   Merges sRNAbench results with blast, featurecounts, miRdeep, miRbase and miRGeneDB.

**Generate all candidates fasta**

path:  /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/

command:

python allCandidatesFasta.py \--all /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/intersections\_table\_sulstoni.xlsx \-s Sulstoni

**Feature Engineering with Ziv’s code**

path: /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Ziv\_Features

command:

python Ziv\_feature\_SOS.py \--precursors /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/all\_candidates\_hairpin.fasta \--mature /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/all\_candidates\_mature.fasta \--star /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/all\_candidates\_star.fasta \--species Sulstoni \--all-remaining /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/intersections\_table\_sulstoni.xlsx

**Statistical Analysis**

Path: /sise/vaksler-group/IsanaRNA/Isana\_Tzah/RNAcentral/miRNAs/Sulstoni/

command:

python statistics.py \--all /sise/vaksler-group/IsanaRNA/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Sulstoni.xlsx \-s Sulstoni