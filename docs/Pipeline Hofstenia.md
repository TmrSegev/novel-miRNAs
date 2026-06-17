**Hofstenia**

**General information:**  
basePath \=  Hofstenia/

Genome and sequencicng libraries were downloaded from:  
https://onedrive.live.com/?authkey=%21AKCytbQPG3KfYmM\&id=746671E5D2B00BDD%2114017\&cid=746671E5D2B00BDD\&sw=bypassConfig

**Preparations:**

1. **Adapter was cut ahead.**  
2. **Combining the trimmed libraries** \- command (done in TrimmedFastq folder):  
   	cat \*fastq \> hofstenia\_final.fastq  
   

**Generating miRNA Candidates Using Mirdeep2**

3. **Genome indexed** by **bowtie-build** \- command:

   bowtie-build \-f Hmia.030120.fasta index/hofsteniaGenomeIndexed

   

   path: \<basePath\>/Hofstenia/genome/index/\*

4. Since we use sequencing data from multiple files we need to **make a config.txt file:**

path: \<basePath\>/bash/config.txt

5. **mapper.pl** \- mapping of the preprocessed reads file to reference database, indexed by bowtie. Split into two parts because of the large size of the libraries. Commands:

	sbatch mapper\_test2.sbatch  
sbatch mapper\_test3.sbatch  
	path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/bash/  
output: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/mapper\_out\_test/	  
	

6. **mirDeep2.pl** \- compared with Hairpin & Mature in all animals.  
   1. Preprocessing command to avoid error of the genome file “has not allowed whitespaces”:  
      	perl \-plane 's/\\s+.+$//' \< caenorhabditis\_Hofstenia.PRJNA13758.WBPS16.genomic.fa \> new\_caenorhabditis\_Hofstenia.PRJNA13758.WBPS16.genomic.fa  
        
   2. **mirDeep2.pl:**  
      Run the mirdeep\_test.sbatch file in each folder separately.  
      path:  
      cd /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/mirdeep\_out/  
      **Or new genome command:**  
      cd /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/mirdeep\_out/  
        
      for dir in AMP1 AMP2 AMP3 DI1 DI2 DI3 EC1 EC2 EC3 GA1 GA2 GA3 HL1 HL2 HL3 IST1 IST2 IST3 PDi1 PDi2 PDi3 PDii1 PDii2 PDii3 PH1 PH2 PH3 PL1 PL2 PL3 SMA1 SMA2 SMA3; do  
        (cd "$dir" && sbatch mirdeep\_test.sbatch)  
      done  
        
        
      

**Generating miRNA Candidates Using sRNABench**

conda activate srnabench

7. **makeSeqObj.jar** \- preprocess genome \- actiavte command:  
   	makeseqobj.sbatch  
   **Full command:**

   java \-jar ../../sRNAtoolboxDB/exec/makeSeqObj.jar /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Genome/refs/Hmia\_ref/Hmia.030120.fasta  
   

   **Command for new genome:**

   java \-jar ../../sRNAtoolboxDB/exec/makeSeqObj.jar /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/sRNA\_PBonly/hofPB\_v6.FINAL.fa  
     
   path: \<basePath\>/sRNAtoolboxDB/seqOBJ/HofsteniaGenomeIndexed.zip  
1. Copying indexed genome to sRNAtoolboxDB/index:

		cp \-r /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Genome/refs/Hmia\_ref/index/. /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB/index/

**Command for new genome:**  
cp \-r /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/sRNA\_PBonly/index/. /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB/index/

2. The seqobj zip file is created in the genome library. Moving command:   
   mv /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Genome/refs/Hmia\_ref/Hmia.zip /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB/seqOBJ/  
   

**Command for new genome:**  
mv /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/sRNA\_PBonly/hofPB\_v6.zip /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB/seqOBJ/

8. **sRNAbench.jar** \- align to genome and prediction for each library \- command:

   java \-jar ../../sRNAtoolboxDB/exec/sRNAbench.jar input=/mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Fastq/Hmia\_annotation/filtered/hofstenia\_final.fastq  output=../../sRNAtoolboxDB/out/Hofstenia predict=true species=hofsteniaGenomeIndexed dbPath=/storage/users/IsanaRNA/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB hairpin=animalsHairpin.fa mature=animalsMature.fa

   

**Filtering**  
	Filter each of the 33 libraries results, trim the sequences and produce a “remaining” file.

**Commands sRNAbench:**  
cd /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/

**new genome:**  
cd /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/scripts/

Filtering and creating remaining.csv files (Note: make sure conda activate is off\!):  
sbatch filter\_hof\_sRNAbench.sbatch

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#  
Instructions: 

1. Run the GFF creation script below once with good\_candidates=False to create the debugging hofstenia file.  
2. Process it using the script process\_debugging\_hofstenia.py.  
3. Then run again the GFF script with good\_candidates=True

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

**Uniting, removing overlaps and creating GFF:**  
python hofsteniasRNAbenchGFF.py \-o Hofstenia\_sRNAbench.gff3 \-seed /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirbase\_data/ALL\_seed\_family\_from\_mirgendb.csv \--create-fasta Hofstenia\_sRNAbench.fasta \-s Hofstenia \--goodcandidates False \--new-genome False

**process\_debugging\_Hofstenia command:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/process\_debugging\_Hofstenia.py \--tool sRNAbench (or miRDeep) (--new-genome)

**New genome:** add boolean flag \--new-genome True

**Commands miRDeep:**  
Filtering and creating remaining.csv files (Note: make sure conda activate is off\!):  
sbatch filter\_hof\_mirdeep.sbatch

python hofsteniaMirdeepFilter.py \-i \*.csv \--filter-s 10 \--filter-mc 10 \--exclude-c 100 

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#  
Instructions: 

4. Run the GFF creation script below once with good\_candidates=False to create the debugging hofstenia file.  
5. Process it using the script process\_debugging\_hofstenia.py.  
6. Then run again the GFF script with good\_candidates=True

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#  
**Uniting, removing overlaps and creating GFF:**  
python hofsteniaMirdeepGFF.py \-o Hofstenia\_mirdeep.gff3 \--create-fasta Hofstenia\_mirdeep.fasta \-s Hofstenia \-seed /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirbase\_data/ALL\_seed\_family\_from\_mirgendb.csv \--goodcandidates False \--new-genome True

**process\_debugging\_Hofstenia command:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/process\_debugging\_Hofstenia.py \--tool miRDeep (--new-genome)

	**New genome:** add boolean flag \--new-genome True

path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/

**Filtering Criteria:**  
There are some pairs of miRNAs that have high scores and have a large overlap \- both should be marked as “overlap”. Marked during the filtering process.

	sRNAbench:

1. **Filter novel.** Remove from novel if:  
   1.  max(5pRC,3pRC)\<10 (weak mature signal)  
   2.  matureBindings\<14 (hairpin doesn’t have enough pairings)  
2. **Filter novel451.**

	**Trimming sequences:**  
	The sequences of the precursor candidates are being trimmed according to the mature/star sequences. That means we cut the beginning and end of the hairpin sequence, so only the parts corresponding to the 5p, loop and 3p sequences are retained. The hairpin sequence and its coordinates are updated in the GFFs. The strand (+ or \-) of the candidate is taken into account.

mirDeep:

1. Discard any sequence that has “rfam alert”.  
2. Discard if a sequence is found in non-coding RNA databases:  
   1. Caenorhabditis\_rRNA.fasta  
   2. Caenorhabditis\_tRNA.fasta  
   3. Caenorhabditis\_snRNA.fasta  
   4. Caenorhabditis\_snoRNA.fasta

   (Location RNAcentral/ncRNAs\_Caenorhabditis/)

3. For sequences with score \<10, check if it is equal to a previous sequence (mature or star) \- if yes \- this is the reason to remove. If score \>10 for both.  
4. Keep a miRNA if it has a score\>=10, OR score\<10 but total\>=100 and star\>0.  
5. Filter if max(mature read count, star read count) \< 10

**For the candidates that are left, we need to mark them as “sense”/”antisense”:**  
“Antisense” miRNAs overlap another miRNA/candidate on the opposite strand:   
The one that has higher counts will be marked as “sense” and the other as “antisense” (the precursor will be marked in the pre-miRNA GFF3 files).  
“Overlap” miRNAs overlap another miRNA/candidate on the **same** strand. Although we filtered overlapping miRNAs earlier, some might be marked as “overlaps” since bedtools calculated overlaps slightly differently than the calculation used for the filter.  
“Sense” and “antisense” strands were found by using bedtools intersect, commands:

1. Preprocessing: remove ‘\\t’ from the end of lines in the gffs.  
   sed \-i 's/\\t\*$//' Hofstenia\_mirdeep\_pre\_only.gff3  
   sed \-i 's/\\t\*$//' Hofstenia\_sRNAbench\_pre\_only.gff3  
     
   path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/  
     
1. Run intersections.sbatch file.  (new genome in Hofstenia\_newGenome)  
   Path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/  
   The commands inside:  
   mirdeep-mirdeep bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_mirdeep\_pre\_only.gff3 \-b /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_mirdeep\_pre\_only.gff3 \> miRdeep\_intersect.bed  
     
   sRNAbench-sRNAbench bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_sRNAbench\_pre\_only.gff3 \-b /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_sRNAbench\_pre\_only.gff3 \> sRNAbench\_intersect.bed  
2. Script commands for marking as overlaps or sense/antisense:  
   mirdeep:  
   python overlapSenseAnti.py \--intersections-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/miRdeep\_intersect.bed \--gff /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_mirdeep\_pre\_only.gff3  
     
   sRNAbench:  
   python overlapSenseAnti.py \--intersections-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/sRNAbench\_intersect.bed \--gff /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_sRNAbench\_pre\_only.gff3  
     
   **New genome:**  
   mirdeep:  
   python overlapSenseAnti.py \--intersections-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/miRdeep\_intersect.bed \--gff /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/scripts/Hofstenia\_mirdeep\_pre\_only.gff3  
     
   sRNAbench:  
   python overlapSenseAnti.py \--intersections-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/sRNAbench\_intersect.bed \--gff /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/scripts/Hofstenia\_sRNAbench\_pre\_only.gff3  
     
   

The output changes the pre\_only gff files in the respective folders.

**Intersecting miRdeep & sRNAbench & Known Results (bedtools)**  
Finding candidates that are part of the intersection between two tools or known miRNAs.  
path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/  
All commands documented in \<path\>/Command.txt

1. When applying intersect force use “-s” for strandness, \-f for minimum overlap.  
   f paramater was chosen based on the “cleanliness” of the canditate. In miRGeneDB and miRDeep the candidates are pures. In sRNAbench, each candidate has an extra head and tail of 11 bases each. miRBase has extra head and tail, sometimes longer.  
   f=0.6 was chosen to allow for the difference between sRNAbench and a pure candidate of roughly 50 bases (50/72=0.694). f=0.5 was chosen to allow for the differences between miRBase and the other databases.  
     
   Apply intersect on:  
   

| sRNAbench | miRDeep | f 0.6 |
| :---- | :---- | :---- |

   

| miRDeep | sRNAbench | f 0.6 |
| :---- | :---- | :---- |

**STAR Alignment**

1)   **Genome indexing** by **STAR** \- command:

   STAR \--runMode genomeGenerate \--runThreadN 16 \--genomeDir ../STAR/genome\_index/ \--genomeSAindexNbases 11  \--genomeFastaFiles /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Genome/refs/Hmia\_ref/Hmia.030120.fasta

   

   command path: \<basePath\>/Hofstenia/bash/star\_genome\_index.sbatch

   output path: \<basePath\>/Hofstenia/STAR/genome\_index/\*

   

2)  **Align** the combined file (including all Hofstenia libraries) **to genome** by **STAR** \- separate command for each library:

   STAR \--genomeDir ../STAR/genome\_index/ \--readFilesIn /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Fastq/Hmia\_annotation/filtered/\<Library\>.filtered.fastq \--outFileNamePrefix ../STAR/align\_to\_genome/\<Library\>/Hofstenia\_ \--outFilterMultimapNmax 20 \--runThreadN 16 \--outSAMtype SAM

	All commands in /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/bash/star\_align\<1 or 2 or 3\>.sbatch

3) **featureCounts** for **sRNAbench gff3 & miRDeep gff3** \- command (run in bash)  
   path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/bash/:  
     
   miRdeep:

featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_mirdeep.gff3 \-o ../counts\_sep/miRNA\_miRdeep\_counts.txt ../STAR/align\_to\_genome/AMP1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA3/Hofstenia\_Aligned.out.sam

	sRNAbench:

featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_sRNAbench.gff3 \-o ../counts\_sep/miRNA\_sRNAbench\_counts.txt ../STAR/align\_to\_genome/AMP1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA3/Hofstenia\_Aligned.out.sam

**For incluster ratio calculation:**

Add flanks of 10bp to pre miRNAs in the GFF:  
python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/add\_flank\_to\_GFF.py (--new-genome)

miRdeep:

featureCounts \-F GFF \-t pre\_miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_mirdeep\_flanked\_pre.gff3 \-o ../counts\_sep/miRNA\_miRdeep\_counts\_flanked.txt ../STAR/align\_to\_genome/AMP1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA3/Hofstenia\_Aligned.out.sam

featureCounts \-F GFF \-t pre\_miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/Hofstenia\_sRNAbench\_flanked\_pre.gff3 \-o ../counts\_sep/miRNA\_sRNAbench\_counts\_flanked.txt ../STAR/align\_to\_genome/AMP1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/AMP3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/DI3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/EC3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/GA3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/HL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/IST3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDi3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PDii3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PH3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/PL3/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA1/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA2/Hofstenia\_Aligned.out.sam ../STAR/align\_to\_genome/SMA3/Hofstenia\_Aligned.out.sam

**Generate Intersections Table**

1. Path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia  
2. Command for script:  
   python intersectionsTableHofstenia.py \-s Hofstenia \--mirdeep-inter-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/miRdeep\_sRNAbench\_intersect.bed \--sRNAbench-inter-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/sRNAbench\_miRdeep\_intersect.bed \--fc-mirdeep /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/counts\_sep/miRNA\_miRdeep\_counts.txt \--fc-pre-mirdeep /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/counts\_sep/miRNA\_miRdeep\_counts\_flanked.txt \--fc-sRNAbench /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/counts\_sep/miRNA\_sRNAbench\_counts.txt \--fc-pre-sRNAbench /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/counts\_sep/miRNA\_sRNAbench\_counts\_flanked.txt \-rm /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/mirdeep\_all\_remaining\_filtered.csv \-rs /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/sRNAbench\_all\_remaining\_filtered.csv \-l EC1,EC2,EC3,GA1,GA2,GA3,DI1,DI2,DI3,PDi1,PDi2,PDi3,PDii1,PDii2,PDii3,PL1,PL2,PL3,PH1,PH2,PH3,HL1,HL2,HL3,IST1,IST2,IST3,AMP1,AMP2,AMP3,SMA1,SMA2,SMA3 \--sum-fc-thres 100  
   Output:  
   intersections\_table\_Hofstenia.xlsx  
     
   **New genome command:**  
   python  /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/intersectionsTableHofstenia.py \-s Hofstenia \--mirdeep-inter-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/miRdeep\_sRNAbench\_intersect.bed \--sRNAbench-inter-table /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/sRNAbench\_miRdeep\_intersect.bed \--fc-mirdeep /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/counts\_sep/miRNA\_miRdeep\_counts.txt \--fc-pre-mirdeep /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/counts\_sep/miRNA\_miRdeep\_counts\_flanked.txt \--fc-sRNAbench /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/counts\_sep/miRNA\_sRNAbench\_counts.txt \--fc-pre-sRNAbench /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/counts\_sep/miRNA\_sRNAbench\_counts\_flanked.txt \-rm /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/scripts/mirdeep\_all\_remaining\_filtered.csv \-rs /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/scripts/sRNAbench\_all\_remaining\_filtered.csv \-l EC1,EC2,EC3,GA1,GA2,GA3,DI1,DI2,DI3,PDi1,PDi2,PDi3,PDii1,PDii2,PDii3,PL1,PL2,PL3,PH1,PH2,PH3,HL1,HL2,HL3,IST1,IST2,IST3,AMP1,AMP2,AMP3,SMA1,SMA2,SMA3 \--sum-fc-thres 100 \--new-genome

**Generate all candidates fasta**

path:  /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/

**command:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/allCandidatesFasta.py \--all /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/intersections\_table\_Hofstenia.xlsx \-s Hofstenia

**New genome command:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/allCandidatesFasta.py \--all /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/intersections\_table\_Hofstenia.xlsx \-s Hofstenia \--new-genome

**Feature Engineering with Ziv’s code**

path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features

**command:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/Ziv\_feature\_SOS.py \--precursors /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/all\_candidates\_hairpin.fasta \--mature /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/all\_candidates\_mature.fasta \--star /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/all\_candidates\_star.fasta \--species Hofstenia \--all-remaining /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/intersections\_table\_Hofstenia.xlsx

**New genome command:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/Ziv\_feature\_SOS.py \--precursors /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/all\_candidates\_hairpin.fasta \--mature /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/all\_candidates\_mature.fasta \--star /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/all\_candidates\_star.fasta \--species Hofstenia\_newGenome \--all-remaining /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia\_newGenome/intersections\_table\_Hofstenia.xlsx \--new-genome

**Feature Engineering Threshold parameters**

Calculating the thresholds for filtering less likely candidates, according to their structural features. Since hofstenia doesn’t have a close relative with known miRNAs, we chose the thresholds according to the known miRNAs of ALL species in mirGenedb. 

Create all fasta files for miRGeneDB

python mirgenedbThresholds.py

path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirbase\_data/

Generate structural features for mirGeneDB

path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/

Run the sbatch file “run.sbatch”.

The command inside the file:

python Ziv\_feature\_SOS.py \--precursors /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirgenedb\_data\_v3/ALL\_mirgenedb\_hairpin.fasta \--mature /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirgenedb\_data\_v3/ALL\_mirgenedb\_mature.fasta \--star /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirgenedb\_data\_v3/ALL\_mirgenedb\_star.fasta \--species miRGeneDB

This script will create temp.csv.

Plot figures

After that run:

python plot\_series.py

to plot the plots from temp.csv

Before choosing the thresholds, we filtered candidates that were regarded as not a valid mir according to the structural features (valid mir \= False), and candidates that had a negative loop. A negative loop length occurred for some candidates with repetitive sequences, in which the mature or the star could occur twice or more in a hairpin sequence. We included calculations of IQR to help determine the thresholds for outlier miRNAs. The lower threshold is calculated as: q25 \- 1.5\*IQR and the upper threshold is calculated as: q75 \+ 1.5\*IQR. (q25=25th percentile, q75=75th percentile)

Hairpin\_seq\_trimmed\_length min: 53 max: 102 lower thres: 55.0 upper thres: 71.0

Mature\_connections min: 14 max: 22 lower thres: 11.5 upper thres: 23.5

Mature\_BP\_ratio min: 0.62 max: 0.96 lower thres: 0.5800000000000001 upper thres: 0.98

Mature\_max\_bulge min: 1.0 max: 4.0 lower thres: \-0.5 upper thres: 3.5

Loop\_length min: 12 max: 54 lower thres: 10.0 upper thres: 26.0

Mature\_Length min: 19 max: 26 lower thres: 20.5 upper thres: 24.5

Star\_length min: 18 max: 25 lower thres: 20.5 upper thres: 24.5

Star\_connections min: 14 max: 23 lower thres: 15.0 upper thres: 23.0

Star\_BP\_ratio min: 0.64 max: 0.96 lower thres: 0.6200000000000001 upper thres: 1.02

Star\_max\_bulge min: 0.0 max: 4.0 lower thres: \-0.5 upper thres: 3.5

Max\_bulge\_symmetry min: 0 max: 3 lower thres: \-1.5 upper thres: 2.5

min\_one\_mer\_hairpin min: 0.1186440677966101 max: 0.2461538461538461 lower thres: 0.1042110874200425 upper thres: 0.27075929874437343

max\_one\_mer\_hairpin min: 0.2542372881355932 max: 0.4098360655737705 lower thres: 0.215507009865737 upper thres: 0.4215154662117053

**Statistical Analysis**

path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/

**command:**

python statistics.py \--all /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Hofstenia.xlsx \-s Hofstenia

The statistics.py script also creates clusters of miRNAs that are within a distance of 10000 bases from each other, on the same strand.

output:

Hofstenia\_clusters\_info.txt \- list of the clusters, with information regarding the seeds and families in each cluster.

hofstenia\_clusters.csv \- only the candidates and their info that are part of a cluster.

**Expression Dynamics**

**For Hofstenia only:**

path: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Expression\_dynamics

**command:** 

python expression\_dynamics.py \--all /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Hofstenia.xlsx \--seed GAGGUAG \--libraries EC1,EC2,EC3,GA1,GA2,GA3,DI1,DI2,DI3,PDi1,PDi2,PDi3,PDii1,PDii2,PDii3,PL1,PL2,PL3,PH1,PH2,PH3,HL1,HL2,HL3,IST1,IST2,IST3,AMP1,AMP2,AMP3,SMA1,SMA2,SMA3 \--time early\_cleavage,gastrula,dimple,pos\_dimple\_phase\_i,post\_dimple\_phase\_ii,pill\_&\_post\_pill,pigmented,pre\_hatchling,hatchling,"in\_situ"\_size\_juvenile,"amputation\_\&RNAi"\_size\_juvenile,sexually\_matured\_adult \-s Hofstenia

**miRge**

The commands below were used either for Oscar’s files, or for the generated files after Ziv’s pipeline, both with old and new genome.

**Prepare for miRge after Ziv:**

**Generate new fasta for remaining candidates after Ziv:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/allCandidatesFasta.py \--all /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Hofstenia.xlsx \-s Hofstenia \--sheetname "(A) Unfiltered" \--output /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/

**New genome:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/RNAcentral/miRNAs/Hofstenia/allCandidatesFasta.py \--all /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Hofstenia\_newGenome.xlsx \-s Hofstenia \--sheetname "(A) Unfiltered" \--output /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/miRge\_after\_Ziv/ 

**Create files necessary for mirge:**

python 	/mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/create\_combined\_mature\_star.py \--base\_path /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/generate\_miRNA\_GFF.py \--species Hofstenia

**new genome:**

python 	/mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/create\_combined\_mature\_star.py \--base\_path /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/miRge\_after\_Ziv/

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/generate\_miRNA\_GFF.py \--species Hofstenia\_newGenome

**(If using Oscar’s files:) Reformat Oscar’s GFF files for mirge-build:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/reformat\_GFF.py \--input /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/miRNA.250203.gff \--output /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/miRNA\_250203\_reformatted.gff3 \--oscar

**Calculating 5p heterogeneity**

Command to run all: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/run\_miRge.sh

Or: /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/run\_miRge.sh

Or one by one (command for after ziv inside run\_miRge.sh in its folder):

conda activate mirge\_env

cd /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_output/

miRge-build   \-g /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Genome/refs/Hmia\_ref/Hmia.030120.fasta   \-mmf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/combined\_mature\_star\_1050.fa   \-hmf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/pre\_miR\_1050\_no\_pre\_in\_seqid.fa   \-mtf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/tRNA.fa   \-ptf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/tRNA.fa   \-snorf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/snoRNA.fa   \-rrf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/rRNA.fa   \-ncof /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/other\_ncRNA.fa   \-mrf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/mRNA.fasta   \-agff /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/miRNA\_no\_pre\_reformatted\_v4.gff3  \-db miRBase   \-on Hofstenia   \-cpu 4

**New genome:**

cd /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/miRge\_after\_Ziv\_output/

miRge-build   \-g /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/sRNA\_PBonly/hofPB\_v6.FINAL.fa   \-mmf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/miRge\_after\_Ziv/combined\_mature\_star.fasta   \-hmf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/miRge\_after\_Ziv/all\_candidates\_hairpin\_T.fasta   \-mtf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/tRNA.fa   \-ptf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/tRNA.fa   \-snorf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/snoRNA.fa   \-rrf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/rRNA.fa   \-ncof /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/other\_ncRNA.fa   \-mrf /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/mRNA.fasta   \-agff /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/miRge\_after\_Ziv/miRNA\_candidates.gff3  \-db miRBase   \-on Hofstenia   \-cpu 4

**Run miRge:**

cd /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/bash

sbatch mirge.sbatch

OR

sbatch mirge\_after\_ziv\_m18.sbatch

**Process miRge results for each library:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/mirge\_processing.py \--dir /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_output/

**After ziv command:**  
python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/mirge\_processing.py \--dir /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv\_output/ \--species Hofstenia \--m18

**New genome:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/mirge\_processing.py \--dir /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia\_newGenome/miRge\_after\_Ziv\_output/ \--species Hofstenia\_newGenome \--m18

**Compare genome to fasta and create good candidates:**

conda activate my\_env

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/compare\_genome\_to\_fasta.py \\  
\--species Hofstenia \\  
\--dir /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge/ \\  
\--genome\_fasta /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Genome/refs/Hmia\_ref/Hmia.030120.fasta \\  
\--gff miRNA\_250203\_reformatted.gff3 \\  
\--premir pre\_miR\_1050\_no\_pre\_in\_seqid.fa \\  
\--mature\_star combined\_mature\_star\_1050.fa \\  
\--output genome\_to\_fasta\_comparison\_new.csv \\  
\--db miRNA\_250203\_reformatted.db

**Command for non-Oscar files but our own generated ones after Ziv:**

python /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/scripts/compare\_genome\_to\_fasta.py \\  
\--species Hofstenia \\  
\--dir /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/miRge\_after\_Ziv/ \\  
\--genome\_fasta /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Hofstenia/Genome/refs/Hmia\_ref/Hmia.030120.fasta \\  
\--gff miRNA\_candidates.gff3 \\  
\--premir all\_candidates\_hairpin\_T.fasta \\  
\--mature\_star combined\_mature\_star.fasta \\  
\--output genome\_to\_fasta\_comparison\_after\_ziv.csv \\  
\--db miRNA\_after\_ziv.db

