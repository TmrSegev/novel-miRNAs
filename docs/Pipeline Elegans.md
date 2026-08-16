**C. Elegans**  

**Scripts directory:** `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/`

The aim of this test is to verify the credibility of the pipeline used on Sultoni and Macrosperma. It is done by using the same pipeline to catalog the MiRNA profile of C. Elegans, and compare it to the already known MiRNA profile.

**General information:**  
basePath \= /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/

Python scripts (novel-miRNAs repo; invoke by absolute path from `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/`):
Per-library filter: srnabenchPerLibraryFilter.py, mirdeepPerLibraryFilter.py  
Unite + GFF + unique_candidates: srnabenchUniteGFF.py, mirdeepUniteGFF.py, processGoodCandidates.py  
Coordinate QC: compare_genome_to_fasta.py (--mode discovery)  
Downstream: overlapSenseAnti.py, add_flank_to_GFF.py, intersectionsTable.py, allCandidatesFasta.py, Ziv_feature_SOS.py, statistics.py  
Legacy combined-run scripts (superseded): sRNAbenchResultsToGFF3.py, mirdeepResultsToGFF3.py

**Preparations:**

1. **Cutting the adapter:**  
   1. Adapter sequence: AACTGTAGGCACCATCAAT  
   2. Sbash file: \<basePath\>/Bash/cutadapt.sbash  
   3. Output summary file: job-2743327.out  
2. **Per-library trimmed FASTQ** (do **not** combine libraries for miRDeep or sRNAbench discovery):  
   Each library is processed separately (CE57 … CE81). Paths: \<basePath\>/TrimmedFastq/SRR\*\_trimmed.fastq
3. **Downloading C. Elegans Genome:**

   wget  [https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis\_elegans/PRJNA13758/caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa.gz](https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS16/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa.gz)

**Text for the paper:**  
The data we used contained small RNA sequencing of three caenorhabditis species: C. Elegans, C. Macrosperma and C. Sulstoni. Small RNA sequencing was downloaded from bioproject PRJNA678899 \[Nelson C, Ambros V. A cohort of Caenorhabditis species lacking the highly conserved let-7 microRNA. G3 (Bethesda). 2021 Apr 23;11(3):jkab022. doi: 10.1093/g3journal/jkab022.  PMID: 33890616; PMCID: PMC8063082.\], and the files were converted from SRA to FASTQ.  
The sequencing libraries:  
C. Sulstoni:

| Run | small RNA-seq of C. sulstoni: L1-arrested L1 larvae on HB101 at 25C |
| :---- | :---- |
| SRR13072575 (SR5) | 12 hours post plating |
| SRR13072574 (SR4) | 16 hours post plating |
| SRR13072573 (SR3) | 20 hours post plating |
| SRR13072572 (SR2) | 24 hours post plating |
| SRR13072571 (SR1) | 28 hours post plating |
| SRR13072570 (SR0) | 32 hours post plating |
| SRR13072577 (SR7) | 4 hours post plating |
| SRR13072576 (SR6) | 8 hours post plating |

C.macrosperma:

| Run | small RNA-seq of C. macrosperma: L1-arrested L1 larvae on HB101 at 25C |
| :---- | :---- |
| SRR13072567 (MR7) | 7 hours post plating |
| SRR13072566 (MR6) | 22 hours post plating |
| SRR13072565 (MR5) | 29 hours post plating |
| SRR13072564 (MR4) | 33 hours post plating |
| SRR13072568 (MR8) | 8 hours post plating |

C.elegans:

| Run | Age | dev\_stage |
| :---- | :---- | :---- |
| SRR13072557 (CE57) | 40 hours post L1 arrest | L3/L4 molt |
| SRR13072558 (CE58) | 36 hours post L1 arrest | L3 |
| SRR13072559 (CE59) | 32 hours post L1 arrest | L3 |
| SRR13072560 (CE60) | 28 hours post L1 arrest | L2/L3 molt |
| SRR13072561 (CE61) | 24 hours post L1 arrest | L2 |
| SRR13072562 (CE62) | 20 hours post L1 arrest | L2 |
| SRR13072563 (CE63) | 16 hours post L1 arrest | L1/L2 molt |
| SRR13072569 (CE69) | 12 hours post L1 arrest | L1 |
| SRR13072578 (CE78) | 48 hours post L1 arrest | L4/adult molt |
| SRR13072579 (CE79) | 44 hours post L1 arrest | L4 |
| SRR13072580 (CE80) | 8 hours post L1 arrest | L1 |
| SRR13072581 (CE81) | 4 hours post L1 arrest | L1 |

The data was processed using cutadapt tool \[DOI:10.14806/ej.17.1.200\]. The adapter sequence (-a AACTGTAGGCACCATCAAT \-e 0.25) was trimmed from each library, while discarding reads that do not contain the adapter, or that are of length smaller than 17 or greater than 26\.

**Generating miRNA Candidates Using Mirdeep2**

4. **Genome indexed** by **bowtie-build** \- command:

   bowtie-build \-f caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa index/elegansGenomeIndexed

   

   path: \<basePath\>/Macrosperma/genome/index/\*

5. **mapper.pl** — Hofstenia-style: **one FASTQ per `mapper.pl` call** (no `config.txt`, no `-d`). Submit from `{base}/Elegans/Bash/`:

   sbatch mapper.sbatch

   Each library writes `mapper_out/elegans_Seq_vs_genome_{LIBRARY}.arf` and `mapper_out/elegans_Seq_collapsed_{LIBRARY}.fasta`.

   Example command (one library):

   mapper.pl .../TrimmedFastq/SRR13072557.1\_trimmed.fastq \-e \-i \-j \-m \-h \-p .../Genome/Index/elegansGenomeIndexed \-t ../mapper\_out/elegans\_Seq\_vs\_genome\_CE57.arf \-s ../mapper\_out/elegans\_Seq\_collapsed\_CE57.fasta

   Flags: `\-e` fastq; `\-i` RNA→DNA; `\-j` drop N; `\-m` collapse; `\-h` fasta out; `\-p` bowtie index; `\-t` ARF; `\-s` collapsed fasta.

6. **mirDeep2.pl** — run **separately for each library** in `{base}/Elegans/mirdeep_out/<library>/` (e.g. `mirdeep_out/CE57/`), using that library’s mapper outputs:

   ```bash
   cd {base}/Elegans/mirdeep_out
   for dir in CE57 CE58 CE59 CE60 CE61 CE62 CE63 CE69 CE78 CE79 CE80 CE81; do
     (cd "$dir" && sbatch mirdeep_test.sbatch)
   done
   # or sequential: sbatch {base}/Elegans/Bash/mirdeep.sbatch
   sbatch {base}/Elegans/scripts/filter_mirdeep.sbatch
   ```

   Outputs per library: `remaining_file_1.csv`, `remaining_file_2.csv`, `removed.csv`.

   **Original combined-run (legacy):** one miRDeep2.pl on unsuffixed `elegans_Seq_collapsed.fasta` into a single `mirdeep_out/` tree — superseded.

7. **mirDeep2.pl** (legacy combined run; see step 6 for per-library workflow) — compared with Hairpin & Mature in all animals.  
   1. Preprocessing command to avoid error of the genome file “has not allowed whitespaces”:  
      	perl \-plane 's/\\s+.+$//' \< caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa \> new\_caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa  
        
   2. mirDeep2.pl command:miRDeep2.pl ../mapper\_out/elegans\_Seq\_collapsed.fasta   
   3. ../Genome/new\_caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa ../mapper\_out/elegans\_Seq\_vs\_genome.arf ../../mirbase\_data/animalsMature.fa none ../../mirbase\_data/animalsHairpin.fa \-g \-1 2\>27.12.2021.mirdeep.log  
      

**Text for the paper:**  
miRNA candidates were generated by using two algorithms: miRDeep2 \[Marc R. Friedländer, Sebastian D. Mackowiak, Na Li, Wei Chen, Nikolaus Rajewsky, miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades, Nucleic Acids Research, Volume 40, Issue 1, 1 January 2012, Pages 37–52, https://doi.org/10.1093/nar/gkr688\] and sRNAbench \[Ernesto Aparicio-Puerta, Ricardo Lebrón, Antonio Rueda, Cristina Gómez-Martín, Stavros Giannoukakos, David Jaspez, José María Medina, Andreja Zubkovic, Igor Jurak, Bastian Fromm, Juan Antonio Marchal, José Oliver, Michael Hackenberg, sRNAbench and sRNAtoolbox 2019: intuitive fast small RNA profiling and differential expression, Nucleic Acids Research, Volume 47, Issue W1, 02 July 2019, Pages W530–W535, https://doi.org/10.1093/nar/gkz415\].  
With the miRDeep algorithm the following steps were taken to generate candidates. The C. Elegans genome \[ref\] was indexed using bowtie-build \[ref\]. Each library FASTQ was mapped with mapper.pl (one command per library; Hofstenia-style). Before the mirDeep2.pl function was used, there was a need for preprocessing to avoid error of the genome file “has not allowed whitespaces”. After that the mirDeep2.pl function was used, compared with hairpin and mature sequences in all animals. Results were generated on 27th Dec. 2021\.

**Generating miRNA Candidates Using sRNABench**

8. **makeSeqObj.jar** \- preprocess genome \- command:

   java \-jar ../../sRNAtoolboxDB/exec/makeSeqObj.jar ../Genome/caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa

   

   path: \<basePath\>/sRNAtoolboxDB/seqOBJ/elegansGenomeIndexed.zip

1. Copying indexed genome to sRNAtoolboxDB/index:

		cp \-r /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/Genome/Index/. /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB/index/

2. The seqobj zip file is created in the genome library. Moving command:   
   mv /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/Genome/caenorhabditis\_elegans.zip /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB/seqOBJ/

9. **sRNAbench.jar** — run **separately for each library** via `sbatch sRNAbench_{LIBRARY}.sbatch` (do not use `elegans_final.fastq`). Example for CE57:

   java \-jar ../../sRNAtoolboxDB/exec/sRNAbench.jar input=../TrimmedFastq/SRR13072557.1\_trimmed.fastq output=../../sRNAtoolboxDB/out/Elegans/Elegans\_CE57 predict=true species=elegansGenomeIndexed dbPath=/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/sRNAtoolboxDB hairpin=animalsHairpin.fa mature=animalsMature.fa

   Repeat for CE58 … CE81 (`output=…/Elegans\_CE58`, etc.). Per-library folders:  
   `{base}/sRNAtoolboxDB/out/Elegans/Elegans_{library}/`

   Filter in each folder (conda off):

   python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py -i novel.txt -a novel451.txt --filter-mc 10

   **Legacy combined run** (superseded): `input=../TrimmedFastq/elegans_final.fastq`, `output=../../sRNAtoolboxDB/out/Elegans`

**Text for the paper:**  
With the sRNAbench algorithm the following steps were taken to generate candidates. First, the genome was preprocessed by making a sequence object. Then the sRNAbench algorithm was used to align to the genome, compared with hairpin and mature sequences in all animals.

**Uniting per-library results, unique_candidates, and creating GFF3/FASTA**

Working directory: `{base}/Elegans/scripts/`  
Per-library filter outputs are read from:  
- sRNAbench: `{base}/sRNAtoolboxDB/out/Elegans/Elegans_{library}/sRNAbench_remaining.csv`  
- miRDeep: `{base}/Elegans/mirdeep_out/{library}/remaining_file_1.csv` and `remaining_file_2.csv`

\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#  
Instructions (same two-pass workflow as Hofstenia):

1. Run each GFF script once with `--uniquecandidates False` → writes `debugging_Elegans_*.csv`  
2. Run `processGoodCandidates.py` for each tool  
3. Run each GFF script again with `--uniquecandidates True` → final GFF + united CSV  
4. Run `compare_genome_to_fasta.py --mode discovery` on each final GFF/FASTA set  
\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#

**sRNAbench — pass 1 (debugging):**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py -o Elegans\_sRNAbench.gff3 -seed /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt --create-fasta Elegans\_sRNAbench.fasta -s Elegans --uniquecandidates False

**processGoodCandidates:**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py --tool sRNAbench -s Elegans

**sRNAbench — pass 2 (final):**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py -o Elegans\_sRNAbench.gff3 -seed /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt --create-fasta Elegans\_sRNAbench.fasta -s Elegans --uniquecandidates True

**miRDeep — pass 1 (debugging):**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py -o Elegans\_mirdeep.gff3 --create-fasta Elegans\_mirdeep.fasta -seed /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt -s Elegans --uniquecandidates False

**processGoodCandidates:**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py --tool miRDeep -s Elegans

**miRDeep — pass 2 (final):**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py -o Elegans\_mirdeep.gff3 --create-fasta Elegans\_mirdeep.fasta -seed /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt -s Elegans --uniquecandidates True

**Coordinate verification (after GFF pass 2):**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/compare_genome_to_fasta.py --mode discovery --species Elegans --dir /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts --genome_fasta /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/Genome/new\_caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa --gff Elegans\_sRNAbench.gff3 --mature Elegans\_sRNAbench.fasta --star Elegans\_sRNAbench\_star.fasta --hairpin-table sRNAbench\_all\_remaining\_filtered.csv --output sRNAbench\_coord\_check.csv

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/compare_genome_to_fasta.py --mode discovery --species Elegans --dir /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts --genome_fasta /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/Genome/new\_caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa --gff Elegans\_mirdeep.gff3 --mature Elegans\_mirdeep.fasta --star Elegans\_mirdeep\_star.fasta --hairpin-table mirdeep\_all\_remaining\_filtered.csv --output mirdeep\_coord\_check.csv

paths (final outputs in scripts/):  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench.gff3  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_mirdeep.gff3  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench\_pre\_only.gff3  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_mirdeep\_pre\_only.gff3  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/sRNAbench\_all\_remaining\_filtered.csv  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/mirdeep\_all\_remaining\_filtered.csv  
/mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/unique\_candidates/

	**Trimming sequences:**  
	Hairpin trimming and coordinate updates are done in **srnabenchPerLibraryFilter.py** (per library, both strands). GFF scripts write trimmed coordinates from the remaining CSV.

	**Filtering Criteria (--filter-mc 10, same as Hofstenia):**

	sRNAbench (srnabenchPerLibraryFilter.py, run per library with `--filter-mc 10`):

1. **Filter novel.** Remove if max(5pRC,3pRC)\<10 or matureBindings\<14.  
2. **Filter novel451.** All novel451 entries are discarded (Hofstenia workflow).  
3. Discard if mature/star/hairpin matches ncRNA database (RNAcentral/ncRNAs\_Caenorhabditis/).  
4. Trim hairpin; drop rows where mature/star cannot be aligned (`sRNAbench_removed_no_find.csv`).

**Unite step (srnabenchUniteGFF.py):** concatenate all libraries → coordinate overlap dedup (≥60%, `overlaps` attribute) → unique_candidates (one representative per ±20 bp cluster; single-library loci kept) → GFF/FASTA.

mirDeep (mirdeepPerLibraryFilter.py, per library: `--filter-s 10 --exclude-c 100 --filter-mc 10`):

1. Discard “rfam alert” or ncRNA match.  
2. Remove lower-scoring duplicate mature/star when score ≥10.  
3. Keep if score≥10, OR score\<10 with total≥100 and star\>0.  
4. Filter if max(mature read count, star read count) \< 10.

**Unite step (mirdeepUniteGFF.py):** same overlap dedup + unique_candidates workflow as sRNAbench.

	MiRDeep Filtering Results (legacy combined run, for reference):

Summary of Input File Number 1 (novel):  
        Total Reads Before Filtering: 636  
        Total Reads Filtered by rfam alert and non-coding RNA database: 15 (2.36%)  
        Total Reads Filtered because they have a duplicate with score \>=10: 95 (15.30%)  
        Total Reads Filtered by Score: 366 (69.58%)  
        Total Reads Filtered because of low mature or star read count: 111 (69.38%)  
        Total Reads Left After all filters: 49 (7.70%)  
Summary of Input File Number 2 (mature):  
        Total Reads Before Filtering: 260  
        Total Reads Filtered by rfam alert and non-coding RNA database: 0 (0.00%)  
        Total Reads Filtered because they have a duplicate with score \>=10: 42 (16.15%)  
        Total Reads Filtered by Score: 35 (16.06%)  
        Total Reads Filtered because of low mature or star read count: 24 (13.11%)  
        Total Reads Left After all filters: 159 (61.15%)

**Text for the paper:**  
sRNAbench candidates were screened per library with `--filter-mc 10`. Sequences were discarded if max(5pRC,3pRC)\<10, matureBindings\<14, or they matched ncRNA databases; all novel451 entries were discarded. Per-library results were united, deduplicated for coordinate overlap, and collapsed to unique_candidates (one representative per ±20 bp cluster; single-library loci kept). miRDeep candidates were filtered per library with the same read-count threshold (10) and score rules as above, then united with the same overlap and unique_candidates collapse.

	**Output files (in Elegans/scripts/):**

1. GFF3 (all features) and `_pre_only.gff3` (precursors only).  
2. Mature, `_pre_only`, and `_star` FASTA files.  
3. `sRNAbench_all_remaining_filtered.csv` / `mirdeep_all_remaining_filtered.csv`.  
4. `debugging_Elegans_*.csv`, `removed_*_no_overlaps.csv`.  
5. Per-library removed CSVs remain in each library folder.

**For the candidates that are left, we need to mark them as “sense”/”antisense” or “overlap”:**  
“Antisense” miRNAs overlap another miRNA/candidate on the **opposite** strand:   
The one that has higher counts will be marked as “sense” and the other as “antisense” (the precursor will be marked in the pre-miRNA GFF3 files).  
“Overlap” miRNAs overlap another miRNA/candidate on the **same** strand.  
Strands were found by using bedtools intersect, commands:  
Path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/

1. Preprocessing: remove ‘\\t’ from the end of lines in the gffs.  
   sed \-i 's/\\t\*$//' Elegans\_mirdeep\_pre\_only.gff3  
   sed \-i 's/\\t\*$//' Elegans\_sRNAbench\_pre\_only.gff3  
2. mirdeep-mirdeep bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_mirdeep\_pre\_only.gff3 \-b /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_mirdeep\_pre\_only.gff3 \> miRdeep\_intersect.bed  
     
   sRNAbench-sRNAbench bedtools intersect command:  
   bedtools intersect \-wao \-loj \-f 0.4 \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench\_pre\_only.gff3 \-b /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench\_pre\_only.gff3 \> sRNAbench\_intersect.bed  
3. Script commands for marking as overlaps or sense/antisense:  
   mirdeep:  
   python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/overlapSenseAnti.py \--intersections-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/miRdeep\_intersect.bed \--gff /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_mirdeep\_pre\_only.gff3  
     
   sRNAbench:  
   python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/overlapSenseAnti.py \--intersections-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench\_intersect.bed \--gff /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench\_pre\_only.gff3  
   

The output changes the pre\_only gff files in the respective folders.  
No overlaps or sense/antisense were found in sRNAbench.

**Text for the paper:**  
Some of the candidates that were left after the filtering were then marked as “sense”/”antisense” or “overlap”. “Overlap” miRNAs overlap another miRNA/candidate on the same strand. “Sense”/“Antisense” miRNAs overlap another miRNA/candidate on the opposite strand.   
The one that has higher counts will be marked as “sense” and the other as “antisense” (the precursor was marked in the pre-miRNA GFF3 files).  
Strands were found by using bedtools intersect \[Aaron R. Quinlan, Ira M. Hall, BEDTools: a flexible suite of utilities for comparing genomic features, Bioinformatics, Volume 26, Issue 6, 15 March 2010, Pages 841–842, https://doi.org/10.1093/bioinformatics/btq033\] (-wao \-loj \-f 0.4) on the same algorithm’s results (for example, sRNAbench results were intersected with themselves).

**Quality Control (mirTrace)**  
[https://github.com/friedlanderlab/mirtrace](https://github.com/friedlanderlab/mirtrace)  
Used as quality control for all libraries.

1. Manual: [https://github.com/friedlanderlab/mirtrace/blob/master/release-bundle-includes/doc/manual/mirtrace\_manual.pdf](https://github.com/friedlanderlab/mirtrace/blob/master/release-bundle-includes/doc/manual/mirtrace_manual.pdf)  
2. Created config file since we use multiple inputs:  
   /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRTrace/  
3. Command for elegans:  
   java \-jar \-Xms4G \-Xmx4G /sise/home/stome/.conda/envs/my\_env/bin/mirtrace.jar qc   \--species cel \--adapter AACTGTAGGCACCATCAAT \--config config.txt

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRTrace/

4. Command for Macrosperma:  
   java \-jar \-Xms4G \-Xmx4G /sise/home/stome/.conda/envs/my\_env/bin/mirtrace.jar qc   \--species cel \--adapter AACTGTAGGCACCATCAAT \--config config.txt

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Macrosperma/miRTrace/

5. Command for Sultoni:  
   java \-jar \-Xms4G \-Xmx4G /sise/home/stome/.conda/envs/my\_env/bin/mirtrace.jar qc   \--species cel \--adapter AACTGTAGGCACCATCAAT \--config config.txt

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Sultoni/miRTrace/

6. Command for Elegans:  
   java \-jar \-Xms4G \-Xmx4G /sise/home/stome/.conda/envs/my\_env/bin/mirtrace.jar qc   \--species cel \--adapter AACTGTAGGCACCATCAAT \--config config.txt

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRTrace/

**Intersecting miRdeep & sRNAbench & mirbase Results (bedtools)**  
Finding candidates that are part of the intersection between two tools or known miRNAs.  
path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans

**Known miRNA GFF files**  
We used Elegans miRNA data taken from two sources: miRbase and miRGeneDB.

miRBase:  
**Editing mirbase GFF**  
Adding the sequences from the fasta files to the GFF file, in addition to trimming the sequences and marking as 5p/3p. Trimming the sequences means we cut the beginning and end of the hairpin sequence, so only the parts corresponding to the 5p, loop and 3p sequences are retained. The hairpin sequence and its coordinates are updated in the GFFs. The strand (+ or \-) of the candidate is taken into account. Candidates with only mature sequence and no star are not trimmed. The mature sequence is only marked as 5p/3p.  
path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/mirbase\_data  
command: python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirbaseToGFF3.py  
output: cel\_mirbase\_seq.gff3

**Text for the paper:**  
We added the sequences from the mirbase fasta files to the mirbase GFF file, in addition to trimming the sequences and marking them as 5p/3p. Trimming the sequences means we cut the beginning and end of the hairpin sequence, so only the parts corresponding to the 5p, loop and 3p sequences are retained. The hairpin sequence and its coordinates were updated in the GFFs. The strand (+ or \-) of the candidate is taken into account. Candidates with only mature sequence and no star are not trimmed. The mature sequence is only marked as 5p/3p.  In the GFFs The “chr” prefix was discarded, and only the number of the chromosome remained.

Download link:  
[https://mirbase.org/ftp/CURRENT/genomes/cel.gff3](https://mirbase.org/ftp/CURRENT/genomes/cel.gff3)  
Preprocessing:   
The “chr” prefix was discarded, and only the number of the chromosome remained.

Mirgendb:  
Download link:  
[https://mirgenedb.org/gff/cel?node=0\&all=1\&sort\_desc=False\&sorted\_by=name\&seed=\&query=\&qtype=mgid\&fnode=0](https://mirgenedb.org/gff/cel?node=0&all=1&sort_desc=False&sorted_by=name&seed=&query=&qtype=mgid&fnode=0)  
Preprocessing:  
There are miRNAs with two versions. Versions 2 were discarded.  
Also, the “chr” prefix was discarded, and only the number of the chromosome remained.

All commands documented in \<path\>/Command.txt  
Run: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/intersections.sbatch

4. When applying intersect force use “-s” for strandness, \-f for minimum overlap.  
   f paramater was chosen based on the “cleanliness” of the canditate. In miRGeneDB and miRDeep the candidates are pures. In sRNAbench, each candidate has an extra head and tail of 11 bases each. miRBase has extra head and tail, sometimes longer.  
   f=0.6 was chosen to allow for the difference between sRNAbench and a pure candidate of roughly 50 bases (50/72=0.694). f=0.5 was chosen to allow for the differences between miRBase and the other databases.  
     
   Apply intersect on:  
   

| sRNAbecnh | miRDeep | f 0.6 |
| :---- | :---- | :---- |
| sRNAbecnh | miRBase | f 0.6 |
| sRNAbecnh | miRGeneDB | f 0.6 |

   

| miRDeep | sRNAbench | f 0.6 |
| :---- | :---- | :---- |
| miRDeep | miRBase | f 0.6 |
| miRDeep | miRGeneDB | f 0.6 |

   

| miRBase | miRGeneDB | f 0.6 |
| :---- | :---- | :---- |
| miRBase | miRDeep | f 0.6 |
| miRBase | sRNAbench | f 0.6 |

**Text for the paper:**  
Next, the results were intersected with each other and with known C. Elegans miRNA databases (miRBase and miRGenedb). It was done in order to increase our certainty regarding the likelihood that a candidate is truly a miRNA, and to test the accuracy of miRdeep and sRNAbench.  
When applying intersect force we used “-s” for strandness, “-f” for minimum overlap.  
“f” parameter was chosen based on the “cleanliness” of the candidate. In miRGeneDB and miRDeep the candidates are pures. In sRNAbench, each candidate has an extra head and tail of 11 bases each. miRBase has extra head and tail, sometimes longer.  
f=0.6 was chosen to allow for the difference between sRNAbench and a pure candidate of roughly 50 bases (50/72=0.694). f=0.5 was chosen to allow for the differences between miRBase and the other databases.  
We created three intersections tables, each focusing on a different set of candidates:

1. Intersecting sRNAbench’s candidates with the others (miRDeep f 0.6, miRBase f 0.5, miRGeneDB f 0.6).  
2. Intersecting miRdeep’s candidates with the others (sRNAbench f 0.6, miRBase f 0.5, miRGeneDB f 0.6).  
3. Intersecting miRBase data with the others  (miRGeneDB f 0.5 miRDeep f 0.5, sRNAbench f 0.5).

   

   

   Results:

   

   **miRBase:**

   Out of 253 miRBase, miRdeep found 155 (61.26%). sRNAbench found 131 (51.77%).

   Out of 138 miRGeneDB, miRdeep found 122 (88.4%). sRNAbench found 109 (78.98%)

   miRdeep seems to be more precise than sRNAbench.

   

   **sRNAbench:**

   Out of 257 results found by sRNAbench (and filtered by us), 130 were known in miRBase and 107 known in miRGeneDB. The last \~100 candidates were not found at all (mostly new mirs).

   58 new mirs were identified. Of those, 8 mirs were also identified by miRDeep,and miRdeep considers 6 of them to be novel. The two “known” are known by miRbase only, not miRGeneDB.

   The 6 novel shared by sRNAbench and miRDeep (as sRNAbench’s name):

   chrV	.	pre\_miRNA	7661068	7661150	.	\-	.	ID=new-mir-novel9;RC\_m=590;RC\_s=1;GACAGGA

   

   chrX	.	pre\_miRNA	14478438	14478538	.	\+	.	ID=new-mir-novel13;RC\_m=168;RC\_s=34;AUCGGUC

   

   chrIV	.	pre\_miRNA	3231126	3231215	.	\-	.	ID=new-mir-novel14;RC\_m=220;RC\_s=3;UCUAGAA

   

   chrIV	.	pre\_miRNA	15400478	15400554	.	\-	.	ID=new-mir-novel15;RC\_m=253;RC\_s=1;AGUCUAA

   

   chrIV	.	pre\_miRNA	11593770	11593860	.	\-	.	ID=new-mir-novel16;RC\_m=129;RC\_s=2;UUCUGGG

   

   chrIV	.	pre\_miRNA	15181491	15181575	.	\+	.	ID=new-mir-novel17;RC\_m=215;RC\_s=1;CUUUUUU

   

   The 2 known shared by sRNAbench and miRDeep (as sRNAbench’s name):

   chrX	.	pre\_miRNA	2329028	2329109	.	\+	.	ID=new-mir-novel6;RC\_m=766;RC\_s=4;UGCCGUA

   

   chrI	.	pre\_miRNA	10600451	10600536	.	\+	.	ID=new-mir-novel12;RC\_m=144;RC\_s=50;cel-mir-2953

	  
The other novels identified by sRNAbench are not identified by miRBase or miRGeneDB.

**miRdeep:**  
Out of 212 results found by miRdeep (and filtered by us), 158 were known in miRBase and 122 known in miRGeneDB.  
49 novel candidates were identified. Out of those, none was known in miRBase or miRGeneDB.  
9 were also identified by sRNAbench, which recognized 3 of them, but 6 are also considered novel by sRNAbench. Those are the same 6 listed above.  
The 3 candidates novel in miRdeep and known in sRNAbench (as miRDeep’s name):  
chrV	.	pre\_miRNA	18903611	18903676	.	\+	.	ID=V\_35152;RC\_m=14071;RC\_s=3;AGUCUUU

chrX	.	pre\_miRNA	1873351	1873408	.	\+	.	ID=X\_83021;RC\_m=612;RC\_s=1;AAAAAGU

chrIII	.	pre\_miRNA	12420399	12420465	.	\-	.	ID=III\_500339;RC\_m=166;RC\_s=27;CAGCGUC

**STAR Alignment**

1)   **Genome indexing** by **STAR** \- command:

   STAR \--runMode genomeGenerate \--runThreadN 16 \--genomeDir ../STAR/genome\_index/ \--genomeSAindexNbases 11  \--genomeFastaFiles ../Genome/caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa

   

   path: \<basePath\>/Elegans/STAR/genome\_index/\*

   

2)  **Align** each library FASTQ **to genome** by **STAR** (per-library; do not use a combined `*_final.fastq`). Command for the first library:

   STAR \--genomeDir ../STAR/genome\_index/ \--readFilesIn ../TrimmedFastq/SRR13072557.1\_trimmed.fastq \--outFileNamePrefix ../STAR/align\_to\_genome/CE57/Elegans\_ \--outFilterMultimapNmax 20 \--runThreadN 16 \--outSAMtype SAM

   

3) **featureCounts** for **sRNAbench gff3 & miRDeep gff3** \- command (run in bash)  
   path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/Bash/:  
     
   miRdeep:

featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a ../scripts/Elegans\_mirdeep.gff3 \-o ../counts\_sep/miRNA\_miRdeep\_counts.txt ../STAR/align\_to\_genome/CE81/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE80/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE69/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE63/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE62/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE61/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE60/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE59/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE58/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE57/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE79/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE78/Elegans\_Aligned.out.sam

	sRNAbench:

featureCounts \-t miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench.gff3 \-o ../counts\_sep/miRNA\_sRNAbench\_counts.txt ../STAR/align\_to\_genome/CE81/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE80/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE69/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE63/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE62/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE61/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE60/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE59/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE58/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE57/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE79/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE78/Elegans\_Aligned.out.sam

mirbase:

featureCounts \-R SAM \-t miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/mirbase\_data/cel\_mirbase\_seq.gff3 \-o ../counts\_sep/miRNA\_mirbase\_counts.txt ../STAR/align\_to\_genome/CE81/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE80/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE69/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE63/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE62/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE61/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE60/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE59/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE58/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE59/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE79/Elegans\_Aligned.out.sam ../STAR/align\_to\_genome/CE78/Elegans\_Aligned.out.sam

**Text for the paper:**  
In order to generate the miRNA counts from our libraries using featurecounts \[ref\], we first had to index the genome using STAR \[ref\] (params?) and align the C. Elegans libraries to the STAR genome. After that, we used featurecounts on each library (-t miRNA \-g ID \-O \-s 1 \-M) for the gff3 files of mirbase, and those generated by miRdeep and sRNAbench.

**For in-cluster m/pre ratio (flanked precursor counts):**

Add flanks of 10 bp to pre_miRNA features in the GFF (run in scripts/):

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/add_flank_to_GFF.py -s Elegans

miRdeep flanked pre_miRNA counts:

featureCounts \-F GFF \-t pre\_miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_mirdeep\_flanked\_pre.gff3 \-o ../counts\_sep/miRNA\_miRdeep\_counts\_flanked.txt \<all STAR SAM files for CE57–CE81\>

sRNAbench flanked pre_miRNA counts:

featureCounts \-F GFF \-t pre\_miRNA \-g ID \-O \-s 1 \-M \-a /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench\_flanked\_pre.gff3 \-o ../counts\_sep/miRNA\_sRNAbench\_counts\_flanked.txt \<all STAR SAM files for CE57–CE81\>

**BLAST** 

1. Remove spaces from the fasta file which will be the blast DB. Command:  
   python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/filterSpacesBlastDB.py \> Caenorhabditis\_pre\_miRNA.fasta  
2. Creating a DB of known miRNAs of nematodes only.  
   1. Create blast DB, command:  
      makeblastdb \-in ../BLAST\_DB/Caenorhabditis\_pre\_miRNA.fasta \-title miRNADB \-dbtype nucl \-out ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB  
3. Blast mature results from miRdeep and sRNAbench. Commands:  
   path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/bash/  
   blastn \-query ../../Charles\_seq/Elegans/scripts/Elegans\_mirdeep.fasta \-db ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB \-out ../queries/Elegans/miRdeep\_blastn\_compact \-outfmt 6 \-evalue 10 \-task blastn-short  
   blastn \-query /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/Elegans\_sRNAbench.fasta \-db ../BLAST\_DB/Caenorhabditis\_pre\_miRNAsDB \-out ../queries/Elegans/sRNAbench\_blastn\_compact \-outfmt 6 \-evalue 10 \-task blastn-short

**Text for the paper:**  
To figure out the closest known homolog miRNA, a BLAST query of the candidates was run. The blast database was created based on the file “Caenorhabditis\_pre\_miRNA.fasta”\[ref\]. Before creating the database, the spaces in the FASTA file were replaced with underscores using the script “filterSpacesBlastDB.py”. The BLAST query’s results were integrated into the intersections tables. Each candidate was matched with its highest scoring matching homolog.

**Generate Intersections Table**

1. Path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/  
2. Command for script:  
   python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/intersectionsTable.py \-s elegans \--mirdeep-inter-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/miRdeep\_sRNAbench\_intersect.bed \--mirdeep-mibrase-inter /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/miRdeep\_miRBase\_intersect.bed \--mirdeep-mirgenedb-inter /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/miRdeep\_miRGeneDB\_intersect.bed \--sRNAbench-inter-table /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench\_miRdeep\_intersect.bed \--sRNAbench-mibrase-inter /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench\_miRBase\_intersect.bed \--sRNAbench-mirgenedb-inter /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench\_miRGeneDB\_intersect.bed \--mirbase-mirgenedb-inter /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/miRBase\_miRGeneDB\_intersect.bed \--mirbase-mirdeep-inter /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/miRBase\_miRdeep\_intersect.bed \--mirbase-sRNAbench-inter /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/miRBase\_sRNAbench\_intersect.bed \--blast-mirdeep /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/queries/Elegans/miRdeep\_blastn\_compact \--blast-sRNAbench /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/queries/Elegans/sRNAbench\_blastn\_compact \--fc-mirdeep /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/counts\_sep/miRNA\_miRdeep\_counts.txt \--fc-pre-mirdeep /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/counts\_sep/miRNA\_miRdeep\_counts\_flanked.txt \--fc-sRNAbench /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/counts\_sep/miRNA\_sRNAbench\_counts.txt \--fc-pre-sRNAbench /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/counts\_sep/miRNA\_sRNAbench\_counts\_flanked.txt \--fc\_mirbase /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/counts\_sep/miRNA\_mirbase\_counts.txt \-rm /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/mirdeep\_all\_remaining\_filtered.csv \-rs /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/sRNAbench\_all\_remaining\_filtered.csv \-mgff /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/mirbase\_data/cel\_mirbase\_seq.gff3 \-l CE81,CE80,CE69,CE63,CE62,CE61,CE60,CE59,CE58,CE57,CE79,CE78 \--sum-fc-thres 100  
   Output:  
   intersections\_table\_script.xlsx, which contains the following (each one in a different sheet):  
* Merges miRDeep results with blast, featurecounts, sRNAbench, miRbase and miRGeneDB.  
* Merges sRNAbench results with blast, featurecounts, miRdeep, miRbase and miRGeneDB.  
* Merges miRBase data  with featurecounts, miRGeneDB, miRdeep and sRNAbench.  
* Blast results for miRdeep.  
* Blast results for sRNAbench.  
3. Every sheet contains categorization by types from 1 to 8\. Each type represents finding by different tools. For example, in mirdeep sheet, type 1 represents candidates found by sRNAbench, miRBase and mirgenedb. Type 2 represents candidates only found by sRNAbench and miRBase.

**Text for the paper:**

For each species, we used a python script to generate Excel files that contain different intersections tables. Each excel file contains the following sheets:

* Intersection of miRDeep as side A, with sRNAbench, miRbase and miRGeneDB as sides B. Also contains BLAST and featurecounts information.  
* Intersection of sRNAbench as side A, with miRDeep, miRbase and miRGeneDB as sides B. Also contains BLAST and featurecounts information.  
* Only for C. Elegans: Intersection of miRBase as side A, with miRGeneDB , miRDeep and sRNAbench as sides B. Also contains BLAST and featurecounts information.  
* Detailed blast results for miRdeep.  
* Detailed blast results for sRNAbench.

Each sheet contains further information regarding the candidates (side A). It contains:

* BLAST information.  
* Featurecounts information.  
* Mature, star and hairpin sequences.  
* Categorization of candidates by types from 1 to 8\. Each type means that the candidate was found by a different combination of tools. For example, in a mirdeep sheet, type 1 represents candidates found by sRNAbench, miRBase and mirgenedb. Type 2 represents candidates only found by sRNAbench and miRBase. Type 8 represents candidates that were only found by miRDeep. Etc.


**Intersections Table Filtering**

Total 253 candidates.

In mirbase sheet, 163/253 have above 100 mature counts in feature counts. 90 were filtered (marked in blue).

Marked in yellow are candidates with a ratio \>2 or \<0.5 of the mature sum featurecounts and the mature readcounts of miRdeep/sRNAbench.

Marked in green are candidates with a ratio \>100 or \<0.01 of the star sum featurecounts and the star readcounts of miRdeep/sRNAbench.

All the novel451 sRNAbench found are considered type 8 (unknown and found only by sRNAbench). Not a single candidate from novel451 intersected with a known candidate. Therefore, **we chose to filter the novel451 from the list of candidates.**

**Intersection Table Manual Changes**

One candidate was changed manually. When intersected as side A it found no intersections, but when intersected as side B it was intersected with both mirdeep and sRNAbench. This candidate is a mirbase candidate that has only mature and no star. Therefore, we did not trim its sequence. The untrimmed sequence was too long for \-f 0.6, and hence the difference between sides A and B.

The mirbase candidate:

ID=MI0010974;Alias=MI0010974;Name=cel-mir-2221;UACAUUUAUUCUUACUGUCCUCGAAUACAAACUGGCGGUUUGCAUUCACUUACAUUUAUAAGACAAAAAUGCAAGUGAUACCAGACCGCUAGUUUGUAAAAGGGAUAAUUUUAUGUGA

Was intersected with miRdeep candidate:

ID=cel-miR-2221;RC\_m=854;RC\_s=44;index=182;MIR-2221

And with sRNAbench candidate:

ID=miR-3781;RC\_m=756;RC\_s=26;index=126;MIR-2221;novel

An with mirgenedb candidate:

ID=Cel-Mir-2221\_pre;Alias=MI0010974

So in the mirbase sheet, we manually added that this mirbase candidate was intersected with the mirdeep, sRNAbench and miRGeneDB candidates, the T/F columns, and also changed its type to type 1\.

In the all candidates sheet, we deleted the mirbase candidate, since it was already part of the mirdeep candidate that intersected with it.

There was another candidate with the same problem, but it belonged to type 5\. It received the same treatment as being updated in the mirbase sheet, and removing the duplicate in the all\_candidates sheet.

The candidate:

ID=MI0000758;Alias=MI0000758;Name=cel-mir-359;AAUGCUCCUUGAAAUUUCAAUCGUUAGAGUAACACACAGUUACACGACCUCAUCAAUCGUGUCACUGGUCUUUCUCUGACGAAUUGAAGUUCUGGAGACAAUUUUGGUUG

In addition, there were 3 candidates that were filtered due to sum\_fc\_m \< 100 in some sheets, but not in others. It caused another disparity. For example, a miRdeep candidate can intersect with a sRNAbench candidate, but this sRNAbench candidate would have been filtered from sRNAbench sheet.  
Those 3 candidates are marked in yellow in the intersections table excel.  
The candidates:  
Type 2:  
ID=MI0010957;Alias=MI0010957;Name=cel-mir-2208b;AAGUGUACCCGGAUCUGAUAUCCUAUCACCAAAAAGAGGAUGCAGAUUUUGGUACACUUCA   
In mirbase sheet, intesected with a mirdeep candidate that was filtered in mirdeep sheet. In type 2 mirbase there are 27 candidates, and in all\_candidates type 2 there are 26 candidates.

ID=MI0008196;Alias=MI0008196;Name=cel-mir-1828;GAUCACUUUUAUCGGUUCCGGUCCCUCUGCAAAAAAGUGGACUGGAAGCAUUUAAGUGAUAGU  
Filtered from mirbase. There are candidates from mirdeep and sRNAbench that are intersected with it.

Type 7:  
ID=miR-4812;RC\_m=4911;RC\_s=6;index=91;AACCACU;novel  
Filtered from sRNAbench. Candidates from mirdeep sheet were intersected with it.

**Text for the paper:**

The intersections tables were filtered based on the featurecounts value. If the sum of featurecounts in all libraries was smaller than 100, the candidate was filtered. This filter corresponds to a decline only for miRdeep candidates, from 211 down to 204\.

Also, we marked outlier candidates in special colors. Marked in yellow are candidates with a ratio \>2 or \<0.5 of the mature sum featurecounts and the mature readcounts of miRdeep/sRNAbench. Marked in green are candidates with a ratio \>100 or \<0.01 of the star sum featurecounts and the star readcounts of miRdeep/sRNAbench.

**Sensitivity Check**

Done on mirbase sheet for C. Elegans as control for our method.

Note: miRGenedb is included in miRBase.

Candidates’ priority:

1. Found by all: miRDeep, sRNAbench, miRBase, miRGenedb (Type 1\)  
2. Found by either miRDeep or sRNAbench and miRBase and miRGenedb (Types 2,3)  
3. In miRGeneDB but not found in miRDeep or sRNAbench. (Type 4\)  
4. Found by miRDeep, sRNAbench and miRBase but not miRGenedb  (Type 5\)  
5. Found by either miRDeep or sRNAbench and miRBase but not miRGenedb (Types 6,7)  
6. In miRBase but not found in miRGenedb, miRDeep or sRNAbench. (Type 8\)  
7. Found by both miRDeep and sRNAbench but is not known by miRBase (novel miRNA) (Type 4 miRDeep / sRNAbench sheets)  
8. Found only by either miRDeep or sRNAbench but not miRBase (Type 8miRDeep / sRNAbench sheets)

**Generate all candidates fasta**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/

command:

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/allCandidatesFasta.py \--all /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/intersections\_table\_elegans.xlsx \-s elegans

**Feature Engineering with Ziv’s code**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Ziv\_Features

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/Ziv_feature_SOS.py \--precursors /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/all\_candidates\_hairpin.fasta \--mature /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/all\_candidates\_mature.fasta \--star /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/all\_candidates\_star.fasta \--species Elegans \--all-remaining /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/intersections\_table\_elegans.xlsx

**Thresholds to filter nematodes**

As detailed in Pipeline Hofstenia, thresholds are chosen from miRGeneDB distributions (run mirgenedbThresholds.py, Ziv with species=miRGeneDB, plot_series.py). Nematodes use the same structural thresholds; sheet (D) also filters **5p_overhang_ziv** and **3p_overhang_ziv** in \[0, 4\]. Description fields use `__` join, `;`→`|`, no periods (set in Ziv_feature_SOS.py and statistics.py).

**Statistical Analysis**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/miRNAs/Elegans/

command:

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/statistics.py \--all /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Elegans.xlsx \-s elegans

**After analyzing all species**

**Interspecies seeds analysis:**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/All\_species

command:

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/seed_frequency.py

**Expression Dynamics**

**For Elegans only:**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/Expression\_dynamics

command: 

run the sbatch file named: **expression\_dynamics.sbatch** (make sure conda is deactivated)

command inside:

xvfb-run \-a python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/expression_dynamics.py \--all /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Elegans.xlsx \--libraries CE81,CE80,CE69,CE63,CE62,CE61,CE60,CE59,CE58,CE57,CE79,CE78 \--time 4,8,12,16,20,24,28,32,36,40,44,48 \-s Elegans

**For all species:**

path: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/All\_species/

run the sbatch file named: **expression\_dynamics\_all.sbatch** (make sure conda is deactivated)

command inside:

xvfb-run \-a python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/expression_dynamics_all.py \--all /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/All\_species/all\_species\_candidates.xlsx

For details and citation purposes see paper "Clustal W and Clustal X

version 2.0", Larkin M., et al. Bioinformatics 2007 23(21):2947-2948

http://bioinformatics.oxfordjournals.org/cgi/content/full/23/21/2947

**miRge**

Calculating 5p heterogeneity.

Generate new fasta for remaining candidates after Ziv:

python /mnt/new_groups/vaksler_group/Isana\_Tzah/novel-miRNAs/allCandidatesFasta.py \--all /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Ziv\_Features/all\_remaining\_after\_ziv\_Elegans.xlsx \-s elegans \--sheetname "(D) Structural Features" \--output /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRge/

Create files necessary for mirge:

python /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/create\_combined\_mature\_star.py

python /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/generate\_miRNA\_GFF.py

If needed, Reformat GFF files for mirge-build:

python /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/reformat\_GFF.py

Command to run all: /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRge/run\_miRge.sh

Or one by one:

conda activate mirge\_env

cd /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRge\_output/

miRge-build   \-g /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/Genome/new\_caenorhabditis\_elegans.PRJNA13758.WBPS16.genomic.fa   \-mmf /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRge/combined\_mature\_star\_1050.fa   \-hmf /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRge/all\_candidates\_hairpin.fasta   \-mtf /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/ncRNAs\_Caenorhabditis/Caenorhabditis\_tRNA.fasta   \-ptf/mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/ncRNAs\_Caenorhabditis/Caenorhabditis\_tRNA.fasta   \-snorf /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/ncRNAs\_Caenorhabditis/Caenorhabditis\_snoRNA.fasta   \-rrf /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/ncRNAs\_Caenorhabditis/Caenorhabditis\_rRNA.fasta   \-ncof /mnt/new_groups/vaksler_group/Isana\_Tzah/RNAcentral/ncRNAs\_Caenorhabditis/Caenorhabditis\_snRNA.fasta   \-mrf /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRge/mRNA.fasta   \-agff /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/miRge/miRNA\_candidates.gff3  \-db miRBase   \-on Elegans   \-cpu 4

cd /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/bash

sbatch mirge.sbatch

Process miRge results for each library:

python /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/mirge\_processing.py

conda activate my\_env

python /mnt/new_groups/vaksler_group/Isana\_Tzah/Charles\_seq/Elegans/scripts/compare\_genome\_to\_fasta.py

---

**New genome track (WBPS19 / PRJNA13758)**

basePath = Elegans\_newGenome/  
Genome: `CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna` under `genome/`  
Annotations: `CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.annotations.gff3` under `genome/`  
Sequencing reads: reuse `Elegans/TrimmedFastq/` (CE57–CE81; same libraries as WBPS16 track)

Download (WBPS19):

wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/caenorhabditis\_elegans/PRJNA13758/caenorhabditis\_elegans.PRJNA13758.WBPS19.genomic.fa.gz

**Preparations (run once from `Elegans_newGenome/bash/`):**

1. **Bowtie index** on WBPS19 scaffolds:

   sbatch bowtie\_index.sbatch

   path: \<basePath\>/genome/index/elegansNewGenomeIndexed.\*

2. **sRNAbench makeSeqObj** on WBPS19:

   sbatch makeseqobj.sbatch

   Move resulting zip to `sRNAtoolboxDB/seqOBJ/elegansNewGenomeIndexed.zip` and copy bowtie index to `sRNAtoolboxDB/index/`.

3. **STAR index** on WBPS19:

   sbatch star\_genome\_indexing.sbatch

   path: \<basePath\>/STAR/genome\_index/

4. **Per-library discovery** (Hofstenia-style; no combined FASTQs):

   sbatch mapper.sbatch  
   sbatch mirdeep.sbatch  
   for lib in CE57 CE58 CE59 CE60 CE61 CE62 CE63 CE69 CE78 CE79 CE80 CE81; do sbatch "sRNAbench_${lib}.sbatch"; done  
   sbatch star\_align.sbatch  
   sbatch ../scripts/filter_mirdeep.sbatch  
   sbatch ../scripts/filter_sRNAbench.sbatch

   Or parallel miRDeep: `for d in CE57 … CE81; do (cd ../mirdeep_out/$d && sbatch mirdeep_test.sbatch); done`

   sRNAbench outputs: `sRNAtoolboxDB/out/Elegans_newGenome/Elegans_{library}/`  
   Do **not** use `star_align_all_libraries.sbatch` (disabled combined legacy).

**Uniting, unique_candidates, and creating GFF3/FASTA (WBPS19 track)**

Working directory: `{base}/Elegans_newGenome/scripts/`  
Use `--variant new_genome` (or `-s Elegans_newGenome` where supported; legacy Hofstenia scripts accept `--new-genome`):

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py -o Elegans\_sRNAbench.gff3 -seed /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt --create-fasta Elegans\_sRNAbench.fasta -s Elegans --variant new\_genome --uniquecandidates False  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py --tool sRNAbench -s Elegans --variant new\_genome  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py -o Elegans\_sRNAbench.gff3 -seed /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt --create-fasta Elegans\_sRNAbench.fasta -s Elegans --variant new\_genome --uniquecandidates True  

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py -o Elegans\_mirdeep.gff3 --create-fasta Elegans\_mirdeep.fasta -seed /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt -s Elegans --variant new\_genome --uniquecandidates False  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py --tool miRDeep -s Elegans --variant new\_genome  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py -o Elegans\_mirdeep.gff3 --create-fasta Elegans\_mirdeep.fasta -seed /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/mirbase\_data/Seeds.txt -s Elegans --variant new\_genome --uniquecandidates True  

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/compare_genome_to_fasta.py --mode discovery --species Elegans --variant new\_genome --dir /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Elegans\_newGenome/scripts --genome\_fasta /mnt/new\_groups/vaksler\_group/Isana\_Tzah/Charles\_seq/Elegans\_newGenome/genome/CELEG.caenorhabditis\_elegans\_PRJNA13758\_WBPS19.scaffolds.fna ...

**Downstream (WBPS19 track)**

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/add_flank_to_GFF.py -s Elegans --variant new\_genome  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/intersectionsTable.py -s Elegans --variant new\_genome ...  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/allCandidatesFasta.py -s Elegans --variant new\_genome --all .../RNAcentral/miRNAs/Elegans\_newGenome/intersections\_table\_Elegans.xlsx  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/Ziv_feature_SOS.py --species Elegans\_newGenome ...  
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/statistics.py -s Elegans\_newGenome --all .../Ziv\_Features/all\_remaining\_after\_ziv\_Elegans\_newGenome.xlsx

Output directory for intersections/FASTAs: `RNAcentral/miRNAs/Elegans_newGenome/`

**miRge (WBPS19 track):** use `Elegans_newGenome/genome/CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna` as `-g` in `miRge-build` instead of the WBPS16 path above.
