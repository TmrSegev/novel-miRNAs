# Novel miRNA discovery pipeline

Unified documentation for all species in this project. Scripts live in `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/`; invoke them by **absolute path** (do not copy into species folders). Data and per-species outputs live under `/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/` and `/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/`.

Species-specific behavior is driven by `pipeline_config.py` (`SPECIES_CONFIG`). Use `-s <Species>` on Python scripts; add `--variant new_genome` (or `-s <Species>_newGenome`) for alternate assemblies.

**Commands in this document are copy-paste ready** with absolute paths. Where a step is only wrapped in sbatch (e.g. per-library STAR align loops), the sbatch file is noted. Species-specific details and paper text remain in `Pipeline Elegans.md`, `Pipeline Hofstenia.md`, `Pipeline Macrosperma.md`, and `Pipeline Sulstoni.md`.

---

## Species at a glance

| Species | Role | Libraries | Discovery mode | Known-miRNA intersects | BLAST | Seed file |
|---------|------|-----------|----------------|----------------------|-------|-----------|
| **Elegans** | Validation control (compare to curated *C. elegans* profile) | 12 (CE57–CE81) | Per-library | **miRBase + miRGeneDB** | Yes | `mirbase_data/Seeds.txt` |
| **Macrosperma** | Novel nematode | 5 (MR4–MR8) | Per-library | No (tool–tool only) | Yes | `mirbase_data/Seeds.txt` |
| **Sulstoni** | Novel nematode | 8 (SR0–SR7) | Per-library | No (tool–tool only) | Yes | `mirbase_data/Seeds.txt` |
| **Hofstenia** | Acoel flatworm | 33 (EC1…SMA3) | Per-library miRDeep; sRNAbench per library | No | **No** | `mirbase_data/ALL_seed_family_from_mirgendb.csv` |

**Nematode sequencing** (Elegans, Macrosperma, Sulstoni): PRJNA678899 (Nelson & Ambros 2021); SRA converted to FASTQ. Adapter: `AACTGTAGGCACCATCAAT`. Trim with cutadapt: `-a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26`.

**Hofstenia**: genome and libraries from [OneDrive dataset](https://onedrive.live.com/?authkey=%21AKCytbQPG3KfYmM&id=746671E5D2B00BDD%2114017&cid=746671E5D2B00BDD&sw=bypassConfig). Adapter trimmed upstream of the pipeline.

### good_candidates support filtering

| Mode | Species | Rule |
|------|---------|------|
| `distinct_libraries` | Elegans, Macrosperma, Sulstoni | ≥2 **distinct libraries** within a 20 bp cluster |
| `dev_condition_replicates` | Hofstenia | ≥2 **developmental-condition replicates** (library name minus trailing digit, e.g. EC1/EC2/EC3 → EC) within a 20 bp cluster |

### Ziv / statistics profiles

| Profile | Species | Primary Ziv input sheet | Extra structural filters |
|---------|---------|-------------------------|--------------------------|
| `structural_sheets` | Nematodes | `(D) Structural Features` | Also filter `5p_overhang_ziv`, `3p_overhang_ziv` ∈ [0, 4] |
| `unfiltered_only` | Hofstenia | `(A) Unfiltered` | Thresholds from pan–miRGeneDB distributions |

Structural thresholds for nematodes follow the same miRGeneDB-derived values as Hofstenia (see [Ziv thresholds](#ziv-structural-thresholds)).

---

## Directory layout

```
Charles_seq/
  <Species>/
    TrimmedFastq/          # nematodes: per-library *_trimmed.fastq
    genome/                # reference FASTA + bowtie index
    mapper_out/            # mapper.pl outputs (nematodes)
    mirdeep_out/<library>/ # per-library miRDeep results
    scripts/               # united GFF/FASTA/CSVs (working dir for unite steps)
    good_candidates/       # processGoodCandidates outputs
    STAR/                  # genome index + per-library SAM
    counts_sep/            # featureCounts outputs
    bash/                  # sbatch wrappers
  <Species>_newGenome/     # alternate assembly track (same library IDs, new indices)
  sRNAtoolboxDB/
    out/<Species>/<Species>_<library>/   # nematode sRNAbench outputs
    out/<Species>_newGenome/...          # new-genome track
  mirbase_data/            # Seeds, miRBase GFF (Elegans), BLAST DB
  Ziv_Features/            # Ziv outputs
RNAcentral/
  miRNAs/<Species>/        # intersections, BED, FASTAs, Excel tables
  queries/<Species>/       # BLAST results (nematodes)
  bash/                    # BLAST batch scripts
```

GFF/FASTA prefixes stay `<Species>_*` even on `{Species}_newGenome` tracks.

---

## Script inventory

| Stage | Script | Notes |
|-------|--------|-------|
| Per-library filter | `srnabenchPerLibraryFilter.py` | `--filter-mc 10` (default) |
| Per-library filter | `mirdeepPerLibraryFilter.py` | `--filter-s 10 --exclude-c 100 --filter-mc 10` |
| Unite + GFF | `srnabenchUniteGFF.py`, `mirdeepUniteGFF.py` | Two-pass `--goodcandidates` |
| good_candidates | `processGoodCandidates.py` | `--tool sRNAbench` or `miRDeep` |
| Coordinate QC | `compare_genome_to_fasta.py` | `--mode discovery`; nematodes (not Hofstenia config default) |
| Strand overlap labels | `overlapSenseAnti.py` | After self-intersect BED |
| Known miRBase GFF | `mirbaseToGFF3.py` | **Elegans only** |
| BLAST DB prep | `filterSpacesBlastDB.py` | Nematodes; once per DB |
| Flanked precursor GFF | `add_flank_to_GFF.py` | 10 bp flanks for m/pre ratio |
| Intersections table | `intersectionsTable.py` | Elegans: miRBase/miRGeneDB args; Hofstenia: no BLAST |
| Candidate FASTAs | `allCandidatesFasta.py` | |
| Structural features | `Ziv_feature_SOS.py` | |
| Statistics / clustering | `statistics.py` | Hofstenia: 10 kb clusters |
| Expression (single species) | `expression_dynamics.py` | Elegans, Hofsteni |
| Expression (all species) | `expression_dynamics_all.py` | After all species complete |
| Seed analysis | `seed_frequency.py` | Cross-species |
| miRGeneDB thresholds | `mirgenedbThresholds.py`, `plot_series.py` | Reference distributions |
| miRge prep | `create_combined_mature_star.py`, `generate_miRNA_GFF.py`, `reformat_GFF.py`, `mirge_processing.py` | Per-species bash wrappers |
| **Legacy (superseded)** | `sRNAbenchResultsToGFF3.py`, `mirdeepResultsToGFF3.py` | Combined-run era |
| **Deprecated wrappers** | `nematode*`, `hofstenia*`, `process_debugging*`, `intersectionsTableHofstenia.py` | Forward to canonical scripts |

---

## Pipeline phases

### 1. Read preprocessing and genome indexing

**Nematodes — cutadapt** (example for one library; repeat per SRR file):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/bash
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 \
  ../Fastq/SRR13072564.1.fastq > ../TrimmedFastq/SRR13072564.1_trimmed.fastq
```

Outputs: `/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/<Species>/TrimmedFastq/*_trimmed.fastq`. Nematodes may also use `bash/cutadapt.sbatch`.

**Nematodes — bowtie index**

```bash
# Elegans
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/genome
bowtie-build -f caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa index/elegansGenomeIndexed

# Macrosperma
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/genome
bowtie-build -f CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna index/macrospermaGenomeIndexed

# Sulstoni
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/genome
bowtie-build -f CSULS.caenorhabditis_sulstoni_JU2788_v1.scaffolds.fna index/sulstoniGenomeIndexed
```

**Nematodes — config.txt** (one trimmed FASTQ per library; example Macrosperma):

```
/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq/SRR13072564.1_trimmed.fastq MR4
/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq/SRR13072565.1_trimmed.fastq MR5
/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq/SRR13072566.1_trimmed.fastq MR6
/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq/SRR13072567.1_trimmed.fastq MR7
/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/TrimmedFastq/SRR13072568.1_trimmed.fastq MR8
```

Path: `/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/<Species>/bash/config.txt`

**Hofstenia**

- Libraries are pre-trimmed upstream. Legacy combined FASTQ (`cat *fastq > hofstenia_final.fastq`) exists but the current workflow filters **per library**.
- Bowtie index:

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/genome
bowtie-build -f Hmia.030120.fasta index/hofsteniaGenomeIndexed
```

- Large mapper runs are split across sbatch jobs (`mapper_test2.sbatch`, `mapper_test3.sbatch`) in `Hofstenia/bash/`.

**Genome whitespace fix** (required before miRDeep2 on some assemblies):

```bash
# Elegans example
perl -lane 's/\s+.+$//' < caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
  > new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa
```

**sRNAbench genome object**

```bash
# Elegans
java -jar /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/exec/makeSeqObj.jar \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/Genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa
cp -r /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/Genome/Index/. \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/index/
mv /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/Genome/caenorhabditis_elegans.zip \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/seqOBJ/elegansGenomeIndexed.zip

# Macrosperma
java -jar /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/exec/makeSeqObj.jar \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna

# Sulstoni
java -jar /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/exec/makeSeqObj.jar \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/genome/CSULS.caenorhabditis_sulstoni_JU2788_v1.scaffolds.fna

# Hofstenia
java -jar /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/exec/makeSeqObj.jar \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/Hmia.030120.fasta
cp -r /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/index/. \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/index/
mv /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/Hmia.zip \
  /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB/seqOBJ/HofsteniaGenomeIndexed.zip
```

---

### 2. miRDeep2 discovery (per library)

**Do not combine FASTQs** for discovery. Each library gets its own folder: `/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/<Species>/mirdeep_out/<library>/`.

**mapper.pl** (nematodes — all libraries via config; run from `<Species>/bash/`):

```bash
# Macrosperma example
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/bash
mapper.pl config.txt -d -e -i -j -m -h \
  -p ../genome/index/macrospermaGenomeIndexed \
  -t ../mapper_out/macrosperma_Seq_vs_genome.arf \
  -s ../mapper_out/macrosperma_Seq_collapsed.fasta

# Sulstoni
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/bash
mapper.pl config.txt -d -e -i -j -m -h \
  -p ../genome/index/sulstoniGenomeIndexed \
  -t ../mapper_out/Sulstoni_Seq_vs_genome.arf \
  -s ../mapper_out/Sulstoni_Seq_collapsed.fasta
```

Flags: `-e` FASTQ; `-d` config; `-p` bowtie prefix; `-i` RNA→DNA; `-j` drop N; `-m` collapse; `-h` FASTA output; `-t` ARF; `-s` collapsed FASTA.

**miRDeep2.pl** — run in each `mirdeep_out/<library>/` (via sbatch or directly), compared against `animalsHairpin.fa` / `animalsMature.fa`.

**Filter** in each library folder (conda off):

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepPerLibraryFilter.py \
  -i result_*.csv --filter-s 10 --exclude-c 100 --filter-mc 10
```

Outputs per library: `remaining_file_1.csv`, `remaining_file_2.csv`, `removed.csv`.

**Hofstenia** — submit miRDeep per library:

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/mirdeep_out
for dir in AMP1 AMP2 AMP3 DI1 DI2 DI3 EC1 EC2 EC3 GA1 GA2 GA3 HL1 HL2 HL3 \
  IST1 IST2 IST3 PDi1 PDi2 PDi3 PDii1 PDii2 PDii3 PH1 PH2 PH3 PL1 PL2 PL3 \
  SMA1 SMA2 SMA3; do
  (cd "$dir" && sbatch mirdeep_test.sbatch)
done
```

Or batch-filter via `sbatch filter_hof_mirdeep.sbatch` in `Hofstenia/scripts/`, then per-library:

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepPerLibraryFilter.py \
  -i *.csv --filter-s 10 --filter-mc 10 --exclude-c 100
```

---

### 3. sRNAbench discovery (per library)

Run **separately for each library** (do not combine FASTQs). Filter in each output folder (conda off).

**Elegans** (example CE57; repeat CE58–CE81):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/bash
java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=../TrimmedFastq/SRR13072557.1_trimmed.fastq \
  output=../../sRNAtoolboxDB/out/Elegans/Elegans_CE57 \
  predict=true species=elegansGenomeIndexed \
  dbPath=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB \
  hairpin=animalsHairpin.fa mature=animalsMature.fa

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py \
  -i novel.txt -a novel451.txt --filter-mc 10
```

**Macrosperma** (example MR4; repeat MR5–MR8):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/bash
java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=../TrimmedFastq/SRR13072564.1_trimmed.fastq \
  output=../../sRNAtoolboxDB/out/Macrosperma/Macrosperma_MR4 \
  predict=true species=macrospermaGenomeIndexed \
  dbPath=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB \
  hairpin=animalsHairpin.fa mature=animalsMature.fa

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py \
  -i novel.txt -a novel451.txt --filter-mc 10
```

**Sulstoni** (example SR0; repeat SR1–SR7):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/bash
java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input=../TrimmedFastq/SRR13072570.1_trimmed.fastq \
  output=../../sRNAtoolboxDB/out/Sulstoni/Sulstoni_SR0 \
  predict=true species=sulstoniGenomeIndexed \
  dbPath=/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/sRNAtoolboxDB \
  hairpin=animalsHairpin.fa mature=animalsMature.fa

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py \
  -i novel.txt -a novel451.txt --filter-mc 10
```

**Hofstenia** — batch filter via `sbatch filter_hof_sRNAbench.sbatch` in `Hofstenia/scripts/`, or per library:

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchPerLibraryFilter.py \
  -i novel.txt -a novel451.txt --filter-mc 10
```

---

### 4. Per-library filtering criteria (`--filter-mc 10`)

**sRNAbench** (`srnabenchPerLibraryFilter.py`)

1. Drop novel if `max(5pRC,3pRC) < 10` or `matureBindings < 14`.
2. **Discard all novel451** entries (all species).
3. Drop if mature/star/hairpin matches ncRNA (`RNAcentral/ncRNAs_Caenorhabditis/`).
4. Trim hairpin to mature/star bounds; drop rows that fail alignment (`sRNAbench_removed_no_find.csv`).

**miRDeep** (`mirdeepPerLibraryFilter.py`)

1. Drop rfam-alert or ncRNA matches.
2. Remove lower-scoring duplicate mature/star when score ≥ 10.
3. Keep if score ≥ 10, **or** score < 10 with total ≥ 100 and star > 0.
4. Drop if `max(mature read count, star read count) < 10`.

---

### 5. Unite libraries, good_candidates, GFF3/FASTA (two-pass)

Working directory: `/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/<Species>/scripts/`

```
1. unite GFF  --goodcandidates False   → debugging CSV
2. processGoodCandidates.py
3. unite GFF  --goodcandidates True    → final GFF + united CSV
4. compare_genome_to_fasta.py --mode discovery   (nematodes)
```

**Elegans — sRNAbench**

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Elegans_sRNAbench.gff3 \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  --create-fasta Elegans_sRNAbench.fasta -s Elegans --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool sRNAbench -s Elegans
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Elegans_sRNAbench.gff3 \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  --create-fasta Elegans_sRNAbench.fasta -s Elegans --goodcandidates True
```

**Elegans — miRDeep**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Elegans_mirdeep.gff3 --create-fasta Elegans_mirdeep.fasta \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  -s Elegans --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool miRDeep -s Elegans
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Elegans_mirdeep.gff3 --create-fasta Elegans_mirdeep.fasta \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  -s Elegans --goodcandidates True
```

**Elegans — coordinate QC**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/compare_genome_to_fasta.py \
  --mode discovery --species Elegans \
  --dir /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts \
  --genome_fasta /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
  --gff Elegans_sRNAbench.gff3 --mature Elegans_sRNAbench.fasta \
  --star Elegans_sRNAbench_star.fasta \
  --hairpin-table sRNAbench_all_remaining_filtered.csv --output sRNAbench_coord_check.csv

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/compare_genome_to_fasta.py \
  --mode discovery --species Elegans \
  --dir /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts \
  --genome_fasta /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/Genome/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa \
  --gff Elegans_mirdeep.gff3 --mature Elegans_mirdeep.fasta \
  --star Elegans_mirdeep_star.fasta \
  --hairpin-table mirdeep_all_remaining_filtered.csv --output mirdeep_coord_check.csv
```

**Macrosperma — sRNAbench and miRDeep** (same two-pass pattern):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Macrosperma_sRNAbench.gff3 \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  --create-fasta Macrosperma_sRNAbench.fasta -s Macrosperma --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool sRNAbench -s Macrosperma
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Macrosperma_sRNAbench.gff3 \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  --create-fasta Macrosperma_sRNAbench.fasta -s Macrosperma --goodcandidates True

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Macrosperma_mirdeep.gff3 --create-fasta Macrosperma_mirdeep.fasta \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  -s Macrosperma --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool miRDeep -s Macrosperma
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Macrosperma_mirdeep.gff3 --create-fasta Macrosperma_mirdeep.fasta \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  -s Macrosperma --goodcandidates True

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/compare_genome_to_fasta.py \
  --mode discovery --species Macrosperma \
  --dir /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts \
  --genome_fasta /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna \
  --gff Macrosperma_sRNAbench.gff3 --mature Macrosperma_sRNAbench.fasta \
  --star Macrosperma_sRNAbench_star.fasta \
  --hairpin-table sRNAbench_all_remaining_filtered.csv --output sRNAbench_coord_check.csv
```

**Sulstoni — sRNAbench and miRDeep**:

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Sulstoni_sRNAbench.gff3 \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  --create-fasta Sulstoni_sRNAbench.fasta -s Sulstoni --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool sRNAbench -s Sulstoni
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Sulstoni_sRNAbench.gff3 \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  --create-fasta Sulstoni_sRNAbench.fasta -s Sulstoni --goodcandidates True

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Sulstoni_mirdeep.gff3 --create-fasta Sulstoni_mirdeep.fasta \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  -s Sulstoni --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool miRDeep -s Sulstoni
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Sulstoni_mirdeep.gff3 --create-fasta Sulstoni_mirdeep.fasta \
  -seed /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Seeds.txt \
  -s Sulstoni --goodcandidates True
```

**Hofstenia — sRNAbench and miRDeep**:

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Hofstenia_sRNAbench.gff3 -s Hofstenia \
  --base-path /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq \
  --create-fasta Hofstenia_sRNAbench.fasta --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool sRNAbench -s Hofstenia \
  --base-path /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/srnabenchUniteGFF.py \
  -o Hofstenia_sRNAbench.gff3 -s Hofstenia \
  --base-path /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq \
  --create-fasta Hofstenia_sRNAbench.fasta --goodcandidates True

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Hofstenia_mirdeep.gff3 --create-fasta Hofstenia_mirdeep.fasta \
  -s Hofstenia --base-path /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq \
  --goodcandidates False
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/processGoodCandidates.py \
  --tool miRDeep -s Hofstenia \
  --base-path /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirdeepUniteGFF.py \
  -o Hofstenia_mirdeep.gff3 --create-fasta Hofstenia_mirdeep.fasta \
  -s Hofstenia --base-path /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq \
  --goodcandidates True
```

**Unite step internals**: concatenate all libraries → coordinate overlap dedup (≥60%, `overlaps` attribute) → good_candidates filter → write GFF/FASTA/`*_all_remaining_filtered.csv`.

**Final outputs** (in `scripts/`):

- `<Species>_sRNAbench.gff3`, `<Species>_mirdeep.gff3` (+ `_pre_only`, `_star`, FASTAs)
- `sRNAbench_all_remaining_filtered.csv`, `mirdeep_all_remaining_filtered.csv`
- `debugging_<Species>_*.csv`, `good_candidates/` artifacts

---

### 6. Sense / antisense / overlap labeling

Precursors that overlap on the genome are labeled using bedtools self-intersection (`-wao -loj -f 0.4 -s`):

| Label | Definition |
|-------|------------|
| **Overlap** | Same-strand overlap ≥ 40% |
| **Sense / antisense** | Opposite-strand overlap; higher-count locus → sense |

**Hofstenia** (run from `Hofstenia/scripts/`; BED output in `RNAcentral/miRNAs/Hofstenia/`):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts
sed -i 's/\t*$//' Hofstenia_mirdeep_pre_only.gff3
sed -i 's/\t*$//' Hofstenia_sRNAbench_pre_only.gff3

cd /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia
bedtools intersect -wao -loj -f 0.4 \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_mirdeep_pre_only.gff3 \
  -b /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_mirdeep_pre_only.gff3 \
  > miRdeep_intersect.bed
bedtools intersect -wao -loj -f 0.4 \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_sRNAbench_pre_only.gff3 \
  -b /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_sRNAbench_pre_only.gff3 \
  > sRNAbench_intersect.bed

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/overlapSenseAnti.py \
  --intersections-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/miRdeep_intersect.bed \
  --gff /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_mirdeep_pre_only.gff3
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/overlapSenseAnti.py \
  --intersections-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/sRNAbench_intersect.bed \
  --gff /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_sRNAbench_pre_only.gff3
```

**Elegans** (same pattern; paths under `Elegans/scripts/` and `RNAcentral/miRNAs/Elegans/`):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts
sed -i 's/\t*$//' Elegans_mirdeep_pre_only.gff3
sed -i 's/\t*$//' Elegans_sRNAbench_pre_only.gff3

cd /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans
bedtools intersect -wao -loj -f 0.4 \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_mirdeep_pre_only.gff3 \
  -b /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_mirdeep_pre_only.gff3 \
  > miRdeep_intersect.bed
bedtools intersect -wao -loj -f 0.4 \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_sRNAbench_pre_only.gff3 \
  -b /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_sRNAbench_pre_only.gff3 \
  > sRNAbench_intersect.bed

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/overlapSenseAnti.py \
  --intersections-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/miRdeep_intersect.bed \
  --gff /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_mirdeep_pre_only.gff3
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/overlapSenseAnti.py \
  --intersections-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench_intersect.bed \
  --gff /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_sRNAbench_pre_only.gff3
```

Macrosperma and Sulstoni follow the same pattern (substitute `Macrosperma` or `Sulstoni` in paths). Or run `intersections.sbatch` in each species folder under `RNAcentral/miRNAs/<Species>/`.

---

### 7. Cross-tool (and known-miRNA) intersections

Use strand-aware intersection (`-s`) with overlap fraction `-f` tuned to precursor “padding”:

| Comparison | `-f` | Rationale |
|------------|------|-----------|
| sRNAbench ↔ miRDeep | 0.6 | sRNAbench adds ~11 bp head/tail (~50/72 ≈ 0.69) |
| Any ↔ miRBase | 0.5–0.6 | miRBase hairpins have variable flanking sequence |
| miRDeep ↔ miRGeneDB | 0.6 | “Pure” coordinates |
| Self (sense/antisense) | 0.4 | Pre-only GFF comparison |

**All nematodes**: sRNAbench ↔ miRDeep only.

**Elegans only** — additional intersections with curated databases:

| Database | Source | Preprocessing |
|----------|--------|---------------|
| **miRBase** | [cel.gff3](https://mirbase.org/ftp/CURRENT/genomes/cel.gff3) | Strip `chr` prefix; build sequence-enriched GFF via `mirbaseToGFF3.py` → `cel_mirbase_seq.gff3` |
| **miRGeneDB** | [cel GFF](https://mirgenedb.org/gff/cel?node=0&all=1) | Drop version-2 duplicates; strip `chr` prefix |

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirbaseToGFF3.py
# output: cel_mirbase_seq.gff3
```

Elegans intersection matrix (all with `-s`):

- sRNAbench ↔ {miRDeep, miRBase, miRGeneDB}
- miRDeep ↔ {sRNAbench, miRBase, miRGeneDB}
- miRBase ↔ {miRGeneDB, miRDeep, sRNAbench}

Full bedtools commands: `RNAcentral/miRNAs/Elegans/Command.txt` and `intersections.sbatch`.

**Elegans validation results** (filtered pipeline vs known DBs):

- miRBase (253): miRDeep 155 (61%); sRNAbench 131 (52%).
- miRGeneDB (138): miRDeep 122 (88%); sRNAbench 109 (79%).
- Shared novel between tools: 6 loci (listed in archived Elegans notes).

---

### 8. STAR alignment and featureCounts

**STAR genome index** (example per species; run from `<Species>/bash/`):

```bash
# Hofstenia
STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir ../STAR/genome_index/ --genomeSAindexNbases 11 \
  --genomeFastaFiles /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/Hmia.030120.fasta

# Elegans
STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir ../STAR/genome_index/ --genomeSAindexNbases 11 \
  --genomeFastaFiles ../Genome/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa

# Macrosperma
STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir ../STAR/genome_index/ --genomeSAindexNbases 11 \
  --genomeFastaFiles ../Genome/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna

# Sulstoni
STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir ../STAR/genome_index/ --genomeSAindexNbases 11 \
  --genomeFastaFiles ../Genome/CSULS.caenorhabditis_sulstoni_JU2788_v1.scaffolds.fna
```

**STAR align** — one command per library (additional libraries in sbatch: `star_align.sbatch`, `star_align1.sbatch`, etc.):

```bash
# Elegans example (CE57)
STAR --genomeDir ../STAR/genome_index/ \
  --readFilesIn ../TrimmedFastq/SRR13072557.1_trimmed.fastq \
  --outFileNamePrefix ../STAR/align_to_genome/CE57/Elegans_ \
  --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM

# Hofstenia example (EC1)
STAR --genomeDir ../STAR/genome_index/ \
  --readFilesIn /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Fastq/Hmia_annotation/filtered/EC1.filtered.fastq \
  --outFileNamePrefix ../STAR/align_to_genome/EC1/Hofstenia_ \
  --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
```

**featureCounts — mature miRNA** (run from `<Species>/bash/`):

**Hofstenia — miRDeep:**

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/bash
featureCounts -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_mirdeep.gff3 \
  -o ../counts_sep/miRNA_miRdeep_counts.txt \
  ../STAR/align_to_genome/AMP1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA3/Hofstenia_Aligned.out.sam
```

**Hofstenia — sRNAbench:**

```bash
featureCounts -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_sRNAbench.gff3 \
  -o ../counts_sep/miRNA_sRNAbench_counts.txt \
  ../STAR/align_to_genome/AMP1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA3/Hofstenia_Aligned.out.sam
```

**Elegans — miRDeep and sRNAbench** (run from `Elegans/bash/`):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/bash
featureCounts -t miRNA -g ID -O -s 1 -M \
  -a ../scripts/Elegans_mirdeep.gff3 \
  -o ../counts_sep/miRNA_miRdeep_counts.txt \
  ../STAR/align_to_genome/CE81/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE80/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE69/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE63/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE62/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE61/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE60/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE59/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE58/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE57/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE79/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE78/Elegans_Aligned.out.sam

featureCounts -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_sRNAbench.gff3 \
  -o ../counts_sep/miRNA_sRNAbench_counts.txt \
  ../STAR/align_to_genome/CE81/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE80/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE69/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE63/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE62/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE61/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE60/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE59/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE58/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE57/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE79/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE78/Elegans_Aligned.out.sam

featureCounts -R SAM -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/cel_mirbase_seq.gff3 \
  -o ../counts_sep/miRNA_mirbase_counts.txt \
  ../STAR/align_to_genome/CE81/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE80/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE69/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE63/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE62/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE61/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE60/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE59/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE58/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE57/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE79/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE78/Elegans_Aligned.out.sam
```

**Macrosperma:**

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/bash
featureCounts -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts/Macrosperma_sRNAbench.gff3 \
  -o ../counts_sep/miRNA_sRNAbench_counts.txt \
  ../STAR/align_to_genome/MR8/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR7/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR6/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR5/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR4/Macrosperma_Aligned.out.sam

featureCounts -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts/Macrosperma_mirdeep.gff3 \
  -o ../counts_sep/miRNA_miRdeep_counts.txt \
  ../STAR/align_to_genome/MR8/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR7/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR6/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR5/Macrosperma_Aligned.out.sam ../STAR/align_to_genome/MR4/Macrosperma_Aligned.out.sam
```

**Sulstoni:**

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/bash
featureCounts -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts/Sulstoni_mirdeep.gff3 \
  -o ../counts_sep/miRNA_miRdeep_counts.txt \
  ../STAR/align_to_genome/SR7/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR6/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR5/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR4/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR3/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR2/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR1/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR0/Sulstoni_Aligned.out.sam

featureCounts -t miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts/Sulstoni_sRNAbench.gff3 \
  -o ../counts_sep/miRNA_sRNAbench_counts.txt \
  ../STAR/align_to_genome/SR7/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR6/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR5/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR4/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR3/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR2/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR1/Sulstoni_Aligned.out.sam ../STAR/align_to_genome/SR0/Sulstoni_Aligned.out.sam
```

**Flanked precursor counts** (m/pre ratio) — add 10 bp flanks, then count `pre_miRNA`:

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/add_flank_to_GFF.py -s <Species>
```

**Hofstenia — flanked miRDeep** (same SAM list as mature counts):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/add_flank_to_GFF.py -s Hofstenia

featureCounts -F GFF -t pre_miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_mirdeep_flanked_pre.gff3 \
  -o ../counts_sep/miRNA_miRdeep_counts_flanked.txt \
  ../STAR/align_to_genome/AMP1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA3/Hofstenia_Aligned.out.sam

featureCounts -F GFF -t pre_miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_sRNAbench_flanked_pre.gff3 \
  -o ../counts_sep/miRNA_sRNAbench_counts_flanked.txt \
  ../STAR/align_to_genome/AMP1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/AMP3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/DI3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/GA3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/HL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/IST3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDi3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PDii3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PH3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/PL3/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA2/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/SMA3/Hofstenia_Aligned.out.sam
```

**Elegans — flanked** (after `add_flank_to_GFF.py -s Elegans`):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/add_flank_to_GFF.py -s Elegans

featureCounts -F GFF -t pre_miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_mirdeep_flanked_pre.gff3 \
  -o ../counts_sep/miRNA_miRdeep_counts_flanked.txt \
  ../STAR/align_to_genome/CE81/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE80/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE69/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE63/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE62/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE61/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE60/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE59/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE58/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE57/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE79/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE78/Elegans_Aligned.out.sam

featureCounts -F GFF -t pre_miRNA -g ID -O -s 1 -M \
  -a /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_sRNAbench_flanked_pre.gff3 \
  -o ../counts_sep/miRNA_sRNAbench_counts_flanked.txt \
  ../STAR/align_to_genome/CE81/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE80/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE69/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE63/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE62/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE61/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE60/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE59/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE58/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE57/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE79/Elegans_Aligned.out.sam ../STAR/align_to_genome/CE78/Elegans_Aligned.out.sam
```

**Macrosperma and Sulstoni** — same flanked pattern after `add_flank_to_GFF.py -s Macrosperma` or `-s Sulstoni`; use the same SAM file lists as the mature-count commands above.

Filter candidates with sum mature featureCounts < 100 (`--sum-fc-thres 100` in intersections table).

---

### 9. BLAST homolog search (nematodes only)

Hofstenia skips BLAST (`use_blast: False`).

**Once — build nematode precursor BLAST DB:**

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/filterSpacesBlastDB.py \
  > Caenorhabditis_pre_miRNA.fasta
makeblastdb -in Caenorhabditis_pre_miRNA.fasta -title miRNADB -dbtype nucl \
  -out /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/Caenorhabditis_pre_miRNAsDB
```

**Per species** (run from `RNAcentral/bash/`):

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/bash

# Elegans
blastn -query ../../Charles_seq/Elegans/scripts/Elegans_mirdeep.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/Elegans/miRdeep_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short
blastn -query /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/Elegans_sRNAbench.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/Elegans/sRNAbench_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short

# Macrosperma
blastn -query ../../Charles_seq/Macrosperma/scripts/Macrosperma_mirdeep.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/Macrosperma/miRdeep_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short
blastn -query /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts/Macrosperma_sRNAbench.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/Macrosperma/sRNAbench_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short

# Sulstoni
blastn -query /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts/Sulstoni_mirdeep.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/Sulstoni/miRdeep_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short
blastn -query /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts/Sulstoni_sRNAbench.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/Sulstoni/sRNAbench_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short
```

---

### 10. Intersections table

**Hofstenia:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/intersectionsTable.py \
  -s Hofstenia \
  --mirdeep-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/miRdeep_sRNAbench_intersect.bed \
  --sRNAbench-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/sRNAbench_miRdeep_intersect.bed \
  --fc-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/counts_sep/miRNA_miRdeep_counts.txt \
  --fc-pre-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/counts_sep/miRNA_miRdeep_counts_flanked.txt \
  --fc-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/counts_sep/miRNA_sRNAbench_counts.txt \
  --fc-pre-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/counts_sep/miRNA_sRNAbench_counts_flanked.txt \
  -rm /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/mirdeep_all_remaining_filtered.csv \
  -rs /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/scripts/sRNAbench_all_remaining_filtered.csv \
  -l EC1,EC2,EC3,GA1,GA2,GA3,DI1,DI2,DI3,PDi1,PDi2,PDi3,PDii1,PDii2,PDii3,PL1,PL2,PL3,PH1,PH2,PH3,HL1,HL2,HL3,IST1,IST2,IST3,AMP1,AMP2,AMP3,SMA1,SMA2,SMA3 \
  --sum-fc-thres 100
```

**Elegans** (includes miRBase, miRGeneDB, BLAST):

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/intersectionsTable.py \
  -s elegans \
  --mirdeep-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/miRdeep_sRNAbench_intersect.bed \
  --mirdeep-mibrase-inter /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/miRdeep_miRBase_intersect.bed \
  --mirdeep-mirgenedb-inter /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/miRdeep_miRGeneDB_intersect.bed \
  --sRNAbench-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench_miRdeep_intersect.bed \
  --sRNAbench-mibrase-inter /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench_miRBase_intersect.bed \
  --sRNAbench-mirgenedb-inter /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/sRNAbench_miRGeneDB_intersect.bed \
  --mirbase-mirgenedb-inter /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/miRBase_miRGeneDB_intersect.bed \
  --mirbase-mirdeep-inter /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/miRBase_miRdeep_intersect.bed \
  --mirbase-sRNAbench-inter /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/miRBase_sRNAbench_intersect.bed \
  --blast-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/queries/Elegans/miRdeep_blastn_compact \
  --blast-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/queries/Elegans/sRNAbench_blastn_compact \
  --fc-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_miRdeep_counts.txt \
  --fc-pre-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_miRdeep_counts_flanked.txt \
  --fc-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_sRNAbench_counts.txt \
  --fc-pre-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_sRNAbench_counts_flanked.txt \
  --fc_mirbase /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_mirbase_counts.txt \
  -rm /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/mirdeep_all_remaining_filtered.csv \
  -rs /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/scripts/sRNAbench_all_remaining_filtered.csv \
  -mgff /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirbase_data/cel_mirbase_seq.gff3 \
  -l CE81,CE80,CE69,CE63,CE62,CE61,CE60,CE59,CE58,CE57,CE79,CE78 \
  --sum-fc-thres 100
```

**Macrosperma:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/intersectionsTable.py \
  -s macrosperma \
  --mirdeep-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Macrosperma/miRdeep_sRNAbench_intersect.bed \
  --sRNAbench-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Macrosperma/sRNAbench_miRdeep_intersect.bed \
  --blast-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/queries/Macrosperma/miRdeep_blastn_compact \
  --blast-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/queries/Macrosperma/sRNAbench_blastn_compact \
  --fc-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/counts_sep/miRNA_miRdeep_counts.txt \
  --fc-pre-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/counts_sep/miRNA_miRdeep_counts_flanked.txt \
  --fc-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/counts_sep/miRNA_sRNAbench_counts.txt \
  --fc-pre-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/counts_sep/miRNA_sRNAbench_counts_flanked.txt \
  -rm /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts/mirdeep_all_remaining_filtered.csv \
  -rs /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Macrosperma/scripts/sRNAbench_all_remaining_filtered.csv \
  -l MR8,MR7,MR6,MR5,MR4 \
  --sum-fc-thres 100
```

**Sulstoni:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/intersectionsTable.py \
  -s sulstoni \
  --mirdeep-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Sulstoni/miRdeep_sRNAbench_intersect.bed \
  --sRNAbench-inter-table /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Sulstoni/sRNAbench_miRdeep_intersect.bed \
  --blast-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/queries/Sulstoni/miRdeep_blastn_compact \
  --blast-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/queries/Sulstoni/sRNAbench_blastn_compact \
  --fc-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/counts_sep/miRNA_miRdeep_counts.txt \
  --fc-pre-mirdeep /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/counts_sep/miRNA_miRdeep_counts_flanked.txt \
  --fc-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/counts_sep/miRNA_sRNAbench_counts.txt \
  --fc-pre-sRNAbench /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/counts_sep/miRNA_sRNAbench_counts_flanked.txt \
  -rm /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts/mirdeep_all_remaining_filtered.csv \
  -rs /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Sulstoni/scripts/sRNAbench_all_remaining_filtered.csv \
  -l SR7,SR6,SR5,SR4,SR3,SR2,SR1,SR0 \
  --sum-fc-thres 100
```

**Output sheets**

| Sheet | Elegans | Other species |
|-------|---------|---------------|
| miRDeep-centric | + miRBase, miRGeneDB, BLAST, FC | + BLAST, FC (nematodes); FC only (Hofstenia) |
| sRNAbench-centric | + miRBase, miRGeneDB, BLAST, FC | same pattern |
| miRBase-centric | **Yes** | No |
| BLAST detail sheets | Nematodes | No |

Candidates are typed 1–8 by which tools/databases recovered them. Elegans supports manual corrections for untrimmed miRBase-only-mature loci and cross-sheet FC filtering edge cases.

---

### 11. Downstream analysis

**Hofstenia — all candidates FASTA:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/allCandidatesFasta.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/intersections_table_Hofstenia.xlsx \
  -s Hofstenia
```

**Hofstenia — Ziv structural features:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/Ziv_feature_SOS.py \
  --precursors /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/all_candidates_hairpin.fasta \
  --mature /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/all_candidates_mature.fasta \
  --star /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/all_candidates_star.fasta \
  --species Hofstenia \
  --all-remaining /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Hofstenia/intersections_table_Hofstenia.xlsx
```

**Hofstenia — statistics:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/statistics.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Hofstenia.xlsx \
  -s Hofstenia
```

**Elegans — all candidates FASTA, Ziv, statistics:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/allCandidatesFasta.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/intersections_table_elegans.xlsx \
  -s elegans

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/Ziv_feature_SOS.py \
  --precursors /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/all_candidates_hairpin.fasta \
  --mature /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/all_candidates_mature.fasta \
  --star /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/all_candidates_star.fasta \
  --species Elegans \
  --all-remaining /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Elegans/intersections_table_elegans.xlsx

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/statistics.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx \
  -s elegans
```

**Macrosperma:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/allCandidatesFasta.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Macrosperma/intersections_table_macrosperma.xlsx \
  -s Macrosperma

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/Ziv_feature_SOS.py \
  --precursors /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Macrosperma/all_candidates_hairpin.fasta \
  --mature /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Macrosperma/all_candidates_mature.fasta \
  --star /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Macrosperma/all_candidates_star.fasta \
  --species Macrosperma \
  --all-remaining /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Macrosperma/intersections_table_macrosperma.xlsx

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/statistics.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Macrosperma.xlsx \
  -s Macrosperma
```

**Sulstoni:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/allCandidatesFasta.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Sulstoni/intersections_table_sulstoni.xlsx \
  -s Sulstoni

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/Ziv_feature_SOS.py \
  --precursors /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Sulstoni/all_candidates_hairpin.fasta \
  --mature /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Sulstoni/all_candidates_mature.fasta \
  --star /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Sulstoni/all_candidates_star.fasta \
  --species Sulstoni \
  --all-remaining /mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/Sulstoni/intersections_table_sulstoni.xlsx

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/statistics.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Sulstoni.xlsx \
  -s Sulstoni
```

Description fields: join with `__`, replace `;` → `|`, strip `ID=` and periods (handled in Ziv/statistics scripts).

---

## Ziv structural thresholds

Derived from all miRGeneDB hairpins. Build reference distributions first:

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/mirgenedbThresholds.py

# Run Ziv on miRGeneDB (command inside Ziv_Features/run.sbatch):
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/Ziv_feature_SOS.py \
  --precursors /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirgenedb_data_v3/ALL_mirgenedb_hairpin.fasta \
  --mature /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirgenedb_data_v3/ALL_mirgenedb_mature.fasta \
  --star /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirgenedb_data_v3/ALL_mirgenedb_star.fasta \
  --species miRGeneDB

python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/plot_series.py
```

Filter invalid miRNAs (`valid_mir = False`) and negative loop lengths before applying IQR bounds (Q25 − 1.5×IQR, Q75 + 1.5×IQR).

| Feature | Lower | Upper |
|---------|-------|-------|
| Hairpin_seq_trimmed_length | 55.0 | 71.0 |
| Mature_connections | 11.5 | 23.5 |
| Mature_BP_ratio | 0.58 | 0.98 |
| Mature_max_bulge | −0.5 | 3.5 |
| Loop_length | 10.0 | 26.0 |
| Mature_Length | 20.5 | 24.5 |
| Star_length | 20.5 | 24.5 |
| Star_connections | 15.0 | 23.0 |
| Star_BP_ratio | 0.62 | 1.02 |
| Star_max_bulge | −0.5 | 3.5 |
| Max_bulge_symmetry | −1.5 | 2.5 |
| min_one_mer_hairpin | 0.104 | 0.271 |
| max_one_mer_hairpin | 0.216 | 0.422 |

Nematodes additionally filter `5p_overhang_ziv` and `3p_overhang_ziv` to [0, 4] on sheet (D).

---

## Species-specific optional steps

### mirTrace QC (nematodes)

[miRTrace](https://github.com/friedlanderlab/mirtrace) on all libraries; `--species cel` proxy, adapter `AACTGTAGGCACCATCAAT`. Config + outputs under `<Species>/miRTrace/`.

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Elegans/miRTrace
java -jar -Xms4G -Xmx4G /sise/home/stome/.conda/envs/my_env/bin/mirtrace.jar qc \
  --species cel --adapter AACTGTAGGCACCATCAAT --config config.txt
```

Same command pattern for Macrosperma and Sulstoni (`Charles_seq/<Species>/miRTrace/`).

### Expression dynamics

**Elegans** (via `expression_dynamics.sbatch`; conda deactivated):

```bash
xvfb-run -a python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/expression_dynamics.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx \
  --libraries CE81,CE80,CE69,CE63,CE62,CE61,CE60,CE59,CE58,CE57,CE79,CE78 \
  --time 4,8,12,16,20,24,28,32,36,40,44,48 \
  -s Elegans
```

**Hofstenia:**

```bash
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/expression_dynamics.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Hofstenia.xlsx \
  --seed GAGGUAG \
  --libraries EC1,EC2,EC3,GA1,GA2,GA3,DI1,DI2,DI3,PDi1,PDi2,PDi3,PDii1,PDii2,PDii3,PL1,PL2,PL3,PH1,PH2,PH3,HL1,HL2,HL3,IST1,IST2,IST3,AMP1,AMP2,AMP3,SMA1,SMA2,SMA3 \
  --time early_cleavage,gastrula,dimple,pos_dimple_phase_i,post_dimple_phase_ii,pill_&_post_pill,pigmented,pre_hatchling,hatchling,"in_situ"_size_juvenile,"amputation_&RNAi"_size_juvenile,sexually_matured_adult \
  -s Hofstenia
```

**All species** (via `expression_dynamics_all.sbatch`):

```bash
xvfb-run -a python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/expression_dynamics_all.py \
  --all /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/All_species/all_species_candidates.xlsx
```

### miRge (5′ heterogeneity)

Run after Ziv filtering. Regenerate FASTAs from the appropriate sheet (`(D) Structural Features` for nematodes, `(A) Unfiltered` for Hofstenia), then `create_combined_mature_star.py`, `generate_miRNA_GFF.py`, `miRge-build`, per-library `mirge.sbatch`, `mirge_processing.py`. Full per-species commands are in the species pipeline docs and `run_miRge.sh` wrappers.

### Elegans sensitivity check

On the miRBase sheet, rank candidates by support across miRDeep, sRNAbench, miRBase, miRGeneDB (types 1–8 priority list in original analysis notes).

### Cross-species seed analysis

After all species complete:

```bash
cd /mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/All_species
python /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/seed_frequency.py
```

---

## Alternate genome assemblies (`--variant new_genome`)

Reuse trimmed reads from the original track; rebuild indices on new scaffolds. Outputs go to `{Species}_newGenome/` with intersections under `RNAcentral/miRNAs/{Species}_newGenome/`.

| Species | New assembly | Genome FASTA |
|---------|--------------|--------------|
| Elegans | WBPS19 / PRJNA13758 | `Elegans_newGenome/genome/CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna` |
| Macrosperma | v2 scaffolds | `Macrosperma_newGenome/genome/CMACR.caenorhabditis_macrosperma_JU2083_v2.scaffolds.fna` |
| Sulstoni | WBPS19 / PRJEB12601 | `Sulstoni_newGenome/genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna` |
| Hofstenia | PacBio | `Hofstenia_newGenome/sRNA_PBonly/hofPB_v6.FINAL.fa` |

**Setup** (from `{Species}_newGenome/bash/`): bowtie index, `makeSeqObj`, STAR index, per-library mapper/miRDeep/sRNAbench/STAR align.

**Python steps**: add `--variant new_genome` to config-aware scripts (or `-s <Species>_newGenome`). Legacy Hofstenia wrappers accept `--new-genome True`.

---

## Library reference (nematodes, PRJNA678899)

### *C. elegans*

| Run | ID | Age (h post L1) | Stage |
|-----|-----|-----------------|-------|
| SRR13072557 | CE57 | 40 | L3/L4 molt |
| SRR13072558 | CE58 | 36 | L3 |
| SRR13072559 | CE59 | 32 | L3 |
| SRR13072560 | CE60 | 28 | L2/L3 molt |
| SRR13072561 | CE61 | 24 | L2 |
| SRR13072562 | CE62 | 20 | L2 |
| SRR13072563 | CE63 | 16 | L1/L2 molt |
| SRR13072569 | CE69 | 12 | L1 |
| SRR13072578 | CE78 | 48 | L4/adult molt |
| SRR13072579 | CE79 | 44 | L4 |
| SRR13072580 | CE80 | 8 | L1 |
| SRR13072581 | CE81 | 4 | L1 |

Genome (original track): WormBase Parasite WBPS16 — `caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa` (whitespace-stripped copy used for miRDeep).

### *C. macrosperma*

| Run | ID | Hours post plating |
|-----|-----|-------------------|
| SRR13072564 | MR4 | 33 |
| SRR13072565 | MR5 | 29 |
| SRR13072566 | MR6 | 22 |
| SRR13072567 | MR7 | 7 |
| SRR13072568 | MR8 | 8 |

Genome: `CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna`

### *C. sulstoni*

| Run | ID | Hours post plating |
|-----|-----|-------------------|
| SRR13072570 | SR0 | 32 |
| SRR13072571 | SR1 | 28 |
| SRR13072572 | SR2 | 24 |
| SRR13072573 | SR3 | 20 |
| SRR13072574 | SR4 | 16 |
| SRR13072575 | SR5 | 12 |
| SRR13072576 | SR6 | 8 |
| SRR13072577 | SR7 | 4 |

Genome: `CSULS.caenorhabditis_sulstoni_JU2788_v1.scaffolds.fna`

### *Hofstenia miamia* (33 libraries)

AMP1–3, DI1–3, EC1–3, GA1–3, HL1–3, IST1–3, PDi1–3, PDii1–3, PH1–3, PL1–3, SMA1–3. Genome: `Hmia.030120.fasta`.

---

## Adding a new species or assembly

1. Add entry to `SPECIES_CONFIG` in `pipeline_config.py` (libraries, seed file, support mode, BLAST/miRBase flags).
2. Place data under `Charles_seq/<Species>/` following the directory layout above.
3. For a new assembly on an existing species, add `NEW_GENOME_OVERRIDES` and use `--variant new_genome`.
4. Run phases 1–11; use Elegans as the validation template if a close relative has curated miRNAs in miRBase/miRGeneDB.

---

## Citations (methods text)

- cutadapt: DOI 10.14806/ej.17.1.200
- miRDeep2: Friedländer et al., NAR 2012 — https://doi.org/10.1093/nar/gkr688
- sRNAbench: Aparicio-Puerta et al., NAR 2019 — https://doi.org/10.1093/nar/gkz415
- bedtools: Quinlan & Hall, Bioinformatics 2010 — https://doi.org/10.1093/bioinformatics/btq033
- Nematode sRNA-seq: Nelson & Ambros, G3 2021 — PRJNA678899
