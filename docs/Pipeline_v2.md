# Novel miRNA discovery pipeline (v2 — template reference)

Generalized, agent-friendly pipeline documentation. Each step gives **one command template** with named placeholders. Substitute values from the [template variables](#template-variables) table and `pipeline_config.py` (`SPECIES_CONFIG`).

| Document | Use when |
|----------|----------|
| **This file** (`Pipeline_v2.md`) | Understanding the workflow; assembling or automating commands |
| [`Pipeline.md`](Pipeline.md) | Copy-paste **full** commands with all library paths expanded |
| `Pipeline <Species>.md` | Species-specific notes, paper text, validation results |

**Paths (fixed):**

- Scripts: `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/`
- Data: `/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/`
- Intersections / BLAST: `/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/`

Invoke Python scripts by **absolute path**; do not copy scripts into species folders.

---

## Template variables

Resolve these before running a step. An agent can load library lists and flags from `pipeline_config.get_species_config("<Species>")`.

| Variable | Meaning | Example (Hofstenia) |
|----------|---------|---------------------|
| `{SPECIES}` | Canonical species name (`-s` argument; **must match** `SPECIES_CONFIG` keys: `Elegans`, `Macrosperma`, `Sulstoni`, `Hofstenia`) | `Hofstenia` |
| `{VARIANT}` | Optional: `--variant new_genome` or omit | *(omit)* |
| `{LIBRARIES}` | Comma-separated library IDs | `EC1,EC2,EC3,...` |
| `{LIBRARY}` | Single library ID | `EC1` |
| `{BASE}` | Charles_seq root | `/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq` |
| `{SPECIES_DIR}` | Species working root | `{BASE}/Hofstenia` |
| `{SCRIPTS_DIR}` | United GFF/FASTA dir | `{SPECIES_DIR}/scripts` |
| `{BASH_DIR}` | sbatch wrappers | `{SPECIES_DIR}/bash` (Elegans: `Bash` with capital B) |
| `{GENOME_DIR}` | Genome folder | `{SPECIES_DIR}/genome` (Elegans: `Genome` with capital G) |
| `{GENOME_FA}` | Reference FASTA | species-specific; see [Library reference](#library-reference) |
| `{GENOME_FA_NO_WS}` | Whitespace-stripped genome (miRDeep) | Elegans: `new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa` |
| `{INDEX_BASENAME}` | bowtie / sRNAbench index basename | `hofstenia`, `elegans`, `macrosperma`, `sulstoni` |
| `{SRNABENCH_INDEX}` | sRNAbench `species=` key | `{INDEX_BASENAME}GenomeIndexed` |
| `{READ_FASTQ}` | Per-library input for STAR/sRNAbench | Nematodes: `{SPECIES_DIR}/TrimmedFastq/<SRR>_trimmed.fastq`; Hofstenia: `{SPECIES_DIR}/Fastq/Hmia_annotation/filtered/{LIBRARY}.filtered.fastq` |
| `{MIRDEEP_OUT}` | Per-library miRDeep folder | `{SPECIES_DIR}/mirdeep_out/{LIBRARY}` |
| `{SRNABENCH_OUT}` | Per-library sRNAbench folder | `{BASE}/sRNAtoolboxDB/out/{SPECIES}/{SPECIES}_{LIBRARY}` |
| `{STAR_INDEX}` | STAR genome dir | `{SPECIES_DIR}/STAR/genome_index` |
| `{STAR_SAM}` | One library SAM | `{SPECIES_DIR}/STAR/align_to_genome/{LIBRARY}/{SPECIES}_Aligned.out.sam` |
| `{STAR_SAMS}` | All library SAMs, space-separated | Expand `{LIBRARY}` over `{LIBRARIES}` |
| `{TOOL}` | Filter/unite script flag | `sRNAbench` or `miRDeep` |
| `{TOOL_TAG}` | Filename / featureCounts token | `sRNAbench` or `mirdeep` (note lowercase **m** in `mirdeep`) |
| `{TOOL_GFF}` | United GFF | `{SCRIPTS_DIR}/{SPECIES}_{TOOL_TAG}.gff3` |
| `{RNACENTRAL}` | RNAcentral root | `/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral` |
| `{SEED}` | Seed file (nematodes only; Hofstenia uses config default) | `{BASE}/mirbase_data/Seeds.txt` |
| `{RNA_MI_DIR}` | Intersections / BED / tables | `/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/{SPECIES}` |
| `{REPO}` | Script root | `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs` |

**Agent assembly rules:**

1. **Per-library loops** — repeat discovery, filter, and STAR-align steps for each ID in `cfg["libraries"]`.
2. **`{STAR_SAMS}`** — join `{STAR_SAM}` for every library:  
   `../STAR/align_to_genome/EC1/Hofstenia_Aligned.out.sam ../STAR/align_to_genome/EC2/...`
3. **`-l` argument** — comma-separated `{LIBRARIES}` with no spaces.
4. **`{VARIANT}`** — append `--variant new_genome` (or `-s {SPECIES}_newGenome`) when using `{Species}_newGenome/` tracks; reuse trimmed reads from the original species folder.
5. **Species forks** — see [Species-specific forks](#species-specific-forks); skip BLAST for Hofstenia; add miRBase/miRGeneDB args for Elegans only.
6. **`-s` on Python scripts** — always use canonical `{SPECIES}` (`Elegans`, not `elegans`). `pipeline_config.get_species_config()` rejects lowercase. Output filenames (e.g. `intersections_table_{SPECIES}.xlsx`) use the `-s` value verbatim.
7. **`{STAR_SAMS}`** — when running featureCounts from `{BASH_DIR}`, prefer relative paths: `../STAR/align_to_genome/{LIBRARY}/{SPECIES}_Aligned.out.sam`.

---

## Species at a glance

| Species | Role | Libraries | Discovery | Known-miRNA intersects | BLAST | Seed file |
|---------|------|-----------|-----------|------------------------|-------|-----------|
| **Elegans** | Validation control | 12 (CE57–CE81) | Per-library | **miRBase + miRGeneDB** | Yes | `mirbase_data/Seeds.txt` |
| **Macrosperma** | Novel nematode | 5 (MR4–MR8) | Per-library | Tool–tool only | Yes | `mirbase_data/Seeds.txt` |
| **Sulstoni** | Novel nematode | 8 (SR0–SR7) | Per-library | Tool–tool only | Yes | `mirbase_data/Seeds.txt` |
| **Hofstenia** | Acoel flatworm | 33 (EC1…SMA3) | Per-library | None | **No** | `mirbase_data/ALL_seed_family_from_mirgendb.csv` |

**Nematodes:** PRJNA678899; adapter `AACTGTAGGCACCATCAAT`; cutadapt  
`-a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26`.

**Hofstenia:** adapter trimmed upstream; large mapper jobs may use sbatch splits.

### good_candidates support

| Mode | Species | Rule |
|------|---------|------|
| `distinct_libraries` | Nematodes | ≥2 distinct libraries in a 20 bp cluster |
| `dev_condition_replicates` | Hofstenia | ≥2 condition replicates (strip trailing digit: EC1/EC2/EC3 → EC) |

### Ziv profiles

| Profile | Species | Input sheet | Extra filters |
|---------|---------|-------------|---------------|
| `structural_sheets` | Nematodes | `(D) Structural Features` | `5p_overhang_ziv`, `3p_overhang_ziv` ∈ [0, 4] |
| `unfiltered_only` | Hofstenia | `(A) Unfiltered` | miRGeneDB-derived thresholds |

---

## Directory layout

```
{BASE}/
  {Species}/
    TrimmedFastq/          # nematodes: *_trimmed.fastq per library
    genome/                # FASTA + bowtie index
    mapper_out/            # mapper.pl (nematodes)
    mirdeep_out/{library}/ # per-library miRDeep
    scripts/               # unite steps, GFF/FASTA/CSVs
    good_candidates/
    STAR/genome_index/
    STAR/align_to_genome/{library}/
    counts_sep/
    bash/                  # sbatch wrappers
  {Species}_newGenome/     # alternate assembly (same library IDs)
  sRNAtoolboxDB/out/{Species}/{Species}_{library}/
  mirbase_data/
  Ziv_Features/
RNAcentral/
  miRNAs/{Species}/        # BED, intersections table, candidate FASTAs
  queries/{Species}/       # BLAST (nematodes)
  bash/
```

GFF/FASTA prefixes stay `{Species}_*` even on `{Species}_newGenome` tracks.

---

## Script inventory

| Stage | Script | Key flags |
|-------|--------|-----------|
| Per-library filter | `srnabenchPerLibraryFilter.py` | `--filter-mc 10` |
| Per-library filter | `mirdeepPerLibraryFilter.py` | `--filter-s 10 --exclude-c 100 --filter-mc 10` |
| Unite + GFF | `srnabenchUniteGFF.py`, `mirdeepUniteGFF.py` | `--goodcandidates False/True` |
| good_candidates | `processGoodCandidates.py` | `--tool sRNAbench` or `miRDeep` |
| Coordinate QC | `compare_genome_to_fasta.py` | `--mode discovery` (nematodes) |
| Strand labels | `overlapSenseAnti.py` | after bedtools self-intersect |
| miRBase GFF | `mirbaseToGFF3.py` | **Elegans only** |
| BLAST DB | `filterSpacesBlastDB.py` | once for nematodes |
| Flanked GFF | `add_flank_to_GFF.py` | `-s {SPECIES}` |
| Intersections table | `intersectionsTable.py` | `--sum-fc-thres 100` |
| Candidate FASTAs | `allCandidatesFasta.py` | |
| Structural features | `Ziv_feature_SOS.py` | |
| Statistics | `statistics.py` | Hofstenia: 10 kb clusters |
| Expression | `expression_dynamics.py`, `expression_dynamics_all.py` | |
| Seed analysis | `seed_frequency.py` | cross-species |

---

## Pipeline phases

### 1. Read preprocessing and genome indexing

**Nematodes — trim adapter** (per library; or `cutadapt.sbatch`):

```bash
cd {BASH_DIR}
cutadapt -a AACTGTAGGCACCATCAAT --core 2 -e 0.25 --discard-untrimmed -m 17 -M 26 \
  ../Fastq/{SRR}.fastq > ../TrimmedFastq/{SRR}_trimmed.fastq
```

**Bowtie index:**

```bash
cd {GENOME_DIR}
bowtie-build -f {GENOME_FA} index/{INDEX_BASENAME}GenomeIndexed
```

**config.txt** (one read file per library for mapper.pl; nematodes only):

```
{READ_FASTQ_1} {LIBRARY_1}
{READ_FASTQ_2} {LIBRARY_2}
...
```

**Genome whitespace fix** (before miRDeep2 if needed):

```bash
perl -lane 's/\s+.+$//' < {GENOME_FA} > {GENOME_FA_NO_WS}
```

**sRNAbench genome object:**

```bash
java -jar {BASE}/sRNAtoolboxDB/exec/makeSeqObj.jar {GENOME_FA}
# Move zip → {BASE}/sRNAtoolboxDB/seqOBJ/{SRNABENCH_INDEX}.zip
# Copy bowtie index → {BASE}/sRNAtoolboxDB/index/
```

> **sbatch:** Hofstenia mapper (`mapper_test2.sbatch`, `mapper_test3.sbatch`); new-genome indexing (`bowtie_index.sbatch`, `makeseqobj.sbatch`, `star_genome_indexing.sbatch`).

---

### 2. miRDeep2 discovery (per library)

**Do not combine FASTQs.** Each library: `{MIRDEEP_OUT}/`.

**mapper.pl** (nematodes — all libraries via config; from `{BASH_DIR}`):

```bash
mapper.pl config.txt -d -e -i -j -m -h \
  -p ../genome/index/{INDEX_BASENAME}GenomeIndexed \
  -t ../mapper_out/{species_tag}_Seq_vs_genome.arf \
  -s ../mapper_out/{species_tag}_Seq_collapsed.fasta
```

**miRDeep2.pl** — run in each `{MIRDEEP_OUT}/` (often via `mirdeep_test.sbatch`).

**Filter** (conda off; inside each library folder):

```bash
python {REPO}/mirdeepPerLibraryFilter.py -i result_*.csv \
  --filter-s 10 --exclude-c 100 --filter-mc 10
```

Outputs: `remaining_file_1.csv`, `remaining_file_2.csv`, `removed.csv`.

**Hofstenia batch submit:**

```bash
cd {SPECIES_DIR}/mirdeep_out
for dir in {LIBRARY_1} {LIBRARY_2} ...; do
  (cd "$dir" && sbatch mirdeep_test.sbatch)
done
```

---

### 3. sRNAbench discovery (per library)

```bash
cd {BASH_DIR}
java -jar ../../sRNAtoolboxDB/exec/sRNAbench.jar \
  input={READ_FASTQ} \
  output=../../sRNAtoolboxDB/out/{SPECIES}/{SPECIES}_{LIBRARY} \
  predict=true species={SRNABENCH_INDEX} \
  dbPath={BASE}/sRNAtoolboxDB \
  hairpin=animalsHairpin.fa mature=animalsMature.fa
```

**Filter** (conda off; in each `{SRNABENCH_OUT}/`):

```bash
python {REPO}/srnabenchPerLibraryFilter.py -i novel.txt -a novel451.txt --filter-mc 10
```

> **sbatch:** `filter_hof_sRNAbench.sbatch`, `filter_hof_mirdeep.sbatch`, `srnabench.sbatch` (new-genome tracks).

---

### 4. Per-library filtering criteria (`--filter-mc 10`)

**sRNAbench:** drop if `max(5pRC,3pRC) < 10` or `matureBindings < 14`; discard all novel451; drop ncRNA matches; trim hairpin to mature/star bounds.

**miRDeep:** drop rfam-alert / ncRNA; deduplicate lower-scoring mature/star when score ≥ 10; keep if score ≥ 10 or (score < 10 and total ≥ 100 and star > 0); drop if `max(mature RC, star RC) < 10`.

---

### 5. Unite libraries, good_candidates, GFF3/FASTA (two-pass)

Working directory: `{SCRIPTS_DIR}/`

```
Pass 1: unite --goodcandidates False  →  debugging CSV
        processGoodCandidates.py
Pass 2: unite --goodcandidates True   →  final GFF + united CSV
        compare_genome_to_fasta.py --mode discovery  (nematodes)
```

**sRNAbench template** (nematodes — pass `-seed`; Hofstenia uses `--base-path` instead):

```bash
cd {SCRIPTS_DIR}

# Pass 1
python {REPO}/srnabenchUniteGFF.py -o {SPECIES}_sRNAbench.gff3 \
  -seed {SEED} --create-fasta {SPECIES}_sRNAbench.fasta \
  -s {SPECIES} {VARIANT} --goodcandidates False

python {REPO}/processGoodCandidates.py --tool sRNAbench -s {SPECIES} {VARIANT}

# Pass 2
python {REPO}/srnabenchUniteGFF.py -o {SPECIES}_sRNAbench.gff3 \
  -seed {SEED} --create-fasta {SPECIES}_sRNAbench.fasta \
  -s {SPECIES} {VARIANT} --goodcandidates True
```

**Hofstenia variant** (no `-seed`; add `--base-path {BASE}`):

```bash
python {REPO}/srnabenchUniteGFF.py -o Hofstenia_sRNAbench.gff3 -s Hofstenia \
  --base-path {BASE} --create-fasta Hofstenia_sRNAbench.fasta --goodcandidates False
python {REPO}/processGoodCandidates.py --tool sRNAbench -s Hofstenia --base-path {BASE}
python {REPO}/srnabenchUniteGFF.py -o Hofstenia_sRNAbench.gff3 -s Hofstenia \
  --base-path {BASE} --create-fasta Hofstenia_sRNAbench.fasta --goodcandidates True
```

**miRDeep** — same two-pass pattern with `mirdeepUniteGFF.py` and `--tool miRDeep`:

```bash
# Nematodes (add -seed; Hofstenia: use --base-path {BASE} instead of -seed)
python {REPO}/mirdeepUniteGFF.py -o {SPECIES}_mirdeep.gff3 \
  --create-fasta {SPECIES}_mirdeep.fasta \
  -seed {SEED} -s {SPECIES} {VARIANT} --goodcandidates False
python {REPO}/processGoodCandidates.py --tool miRDeep -s {SPECIES} {VARIANT}
python {REPO}/mirdeepUniteGFF.py -o {SPECIES}_mirdeep.gff3 \
  --create-fasta {SPECIES}_mirdeep.fasta \
  -seed {SEED} -s {SPECIES} {VARIANT} --goodcandidates True
```

**Hofstenia miRDeep unite** (add `--base-path {BASE}` on all three steps; omit `-seed`):

```bash
python {REPO}/mirdeepUniteGFF.py -o Hofstenia_mirdeep.gff3 \
  --create-fasta Hofstenia_mirdeep.fasta -s Hofstenia \
  --base-path {BASE} --goodcandidates False
python {REPO}/processGoodCandidates.py --tool miRDeep -s Hofstenia --base-path {BASE}
python {REPO}/mirdeepUniteGFF.py -o Hofstenia_mirdeep.gff3 \
  --create-fasta Hofstenia_mirdeep.fasta -s Hofstenia \
  --base-path {BASE} --goodcandidates True
```

**Coordinate QC** (nematodes; use `{GENOME_FA_NO_WS}` for Elegans):

```bash
python {REPO}/compare_genome_to_fasta.py --mode discovery --species {SPECIES} {VARIANT} \
  --dir {SCRIPTS_DIR} --genome_fasta {GENOME_FA_NO_WS} \
  --gff {SPECIES}_{TOOL_TAG}.gff3 --mature {SPECIES}_{TOOL_TAG}.fasta \
  --star {SPECIES}_{TOOL_TAG}_star.fasta \
  --hairpin-table {TOOL_TAG}_all_remaining_filtered.csv --output {TOOL_TAG}_coord_check.csv
```

**Final outputs in `{SCRIPTS_DIR}/`:** `{SPECIES}_sRNAbench.gff3`, `{SPECIES}_mirdeep.gff3`, `*_pre_only.gff3`, FASTAs, `*_all_remaining_filtered.csv`, `debugging_{SPECIES}_*.csv`.

---

### 6. Sense / antisense / overlap labeling

Self-intersection on `_pre_only.gff3` files (`bedtools intersect -wao -loj -f 0.4`):

```bash
cd {SCRIPTS_DIR}
sed -i 's/\t*$//' {SPECIES}_mirdeep_pre_only.gff3
sed -i 's/\t*$//' {SPECIES}_sRNAbench_pre_only.gff3

cd {RNA_MI_DIR}
bedtools intersect -wao -loj -f 0.4 \
  -a {SCRIPTS_DIR}/{SPECIES}_mirdeep_pre_only.gff3 \
  -b {SCRIPTS_DIR}/{SPECIES}_mirdeep_pre_only.gff3 \
  > miRdeep_intersect.bed
bedtools intersect -wao -loj -f 0.4 \
  -a {SCRIPTS_DIR}/{SPECIES}_sRNAbench_pre_only.gff3 \
  -b {SCRIPTS_DIR}/{SPECIES}_sRNAbench_pre_only.gff3 \
  > sRNAbench_intersect.bed

python {REPO}/overlapSenseAnti.py \
  --intersections-table {RNA_MI_DIR}/miRdeep_intersect.bed \
  --gff {SCRIPTS_DIR}/{SPECIES}_mirdeep_pre_only.gff3
python {REPO}/overlapSenseAnti.py \
  --intersections-table {RNA_MI_DIR}/sRNAbench_intersect.bed \
  --gff {SCRIPTS_DIR}/{SPECIES}_sRNAbench_pre_only.gff3
```

> **sbatch:** `intersections.sbatch` in `{RNA_MI_DIR}/` runs bedtools for cross-tool intersections.

---

### 7. Cross-tool (and known-miRNA) intersections

Strand-aware bedtools (`-s`); overlap fraction `-f`:

| Comparison | `-f` |
|------------|------|
| sRNAbench ↔ miRDeep | 0.6 |
| Any ↔ miRBase | 0.5–0.6 |
| miRDeep ↔ miRGeneDB | 0.6 |
| Self (sense/antisense) | 0.4 |

**Nematodes:** sRNAbench ↔ miRDeep only.

**Elegans only** — build miRBase GFF, then intersect with miRBase + miRGeneDB:

```bash
cd {BASE}/mirbase_data
python {REPO}/mirbaseToGFF3.py   # → cel_mirbase_seq.gff3
```

Full bedtools commands: `{RNA_MI_DIR}/Command.txt` and `intersections.sbatch`.

---

### 8. STAR alignment and featureCounts

**STAR index** (from `{BASH_DIR}`):

```bash
STAR --runMode genomeGenerate --runThreadN 16 \
  --genomeDir ../STAR/genome_index/ --genomeSAindexNbases 11 \
  --genomeFastaFiles {GENOME_FA}
```

**STAR align** (one command per `{LIBRARY}`; remaining libraries in `star_align*.sbatch`):

```bash
STAR --genomeDir ../STAR/genome_index/ \
  --readFilesIn {READ_FASTQ} \
  --outFileNamePrefix ../STAR/align_to_genome/{LIBRARY}/{SPECIES}_ \
  --outFilterMultimapNmax 20 --runThreadN 16 --outSAMtype SAM
```

**featureCounts — mature miRNA** (from `{BASH_DIR}`; expand `{STAR_SAMS}` over all libraries):

```bash
cd {BASH_DIR}
featureCounts -t miRNA -g ID -O -s 1 -M \
  -a {SCRIPTS_DIR}/{SPECIES}_{TOOL_TAG}.gff3 \
  -o ../counts_sep/miRNA_{TOOL_TAG}_counts.txt \
  {STAR_SAMS}
```

**Elegans only — miRBase counts:**

```bash
featureCounts -R SAM -t miRNA -g ID -O -s 1 -M \
  -a {BASE}/mirbase_data/cel_mirbase_seq.gff3 \
  -o ../counts_sep/miRNA_mirbase_counts.txt \
  {STAR_SAMS}
```

**Flanked precursor counts** (m/pre ratio):

```bash
python {REPO}/add_flank_to_GFF.py -s {SPECIES} {VARIANT}

featureCounts -F GFF -t pre_miRNA -g ID -O -s 1 -M \
  -a {SCRIPTS_DIR}/{SPECIES}_{TOOL_TAG}_flanked_pre.gff3 \
  -o ../counts_sep/miRNA_{TOOL_TAG}_counts_flanked.txt \
  {STAR_SAMS}
```

Filter candidates with sum mature featureCounts < 100 (`--sum-fc-thres 100` in intersections table).

> **Worked example:** see [`Pipeline.md` §8](Pipeline.md) or `Pipeline Hofstenia.md` for Hofstenia with all 33 SAM paths expanded.

---

### 9. BLAST homolog search (nematodes only)

Hofstenia: `use_blast: False` — skip this phase.

**Once — build DB** (output lives in `{RNACENTRAL}/BLAST_DB/`, matching the `blastn -db` path):

```bash
cd {RNACENTRAL}/BLAST_DB
python {REPO}/filterSpacesBlastDB.py > Caenorhabditis_pre_miRNA.fasta
makeblastdb -in Caenorhabditis_pre_miRNA.fasta -title miRNADB -dbtype nucl \
  -out Caenorhabditis_pre_miRNAsDB
```

**Per species** (from `{RNACENTRAL}/bash/`):

```bash
blastn -query {SCRIPTS_DIR}/{SPECIES}_mirdeep.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/{SPECIES}/miRdeep_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short
blastn -query {SCRIPTS_DIR}/{SPECIES}_sRNAbench.fasta \
  -db ../BLAST_DB/Caenorhabditis_pre_miRNAsDB \
  -out ../queries/{SPECIES}/sRNAbench_blastn_compact \
  -outfmt 6 -evalue 10 -task blastn-short
```

---

### 10. Intersections table

**Base template** (all species):

```bash
python {REPO}/intersectionsTable.py -s {SPECIES} {VARIANT} \
  --mirdeep-inter-table {RNA_MI_DIR}/miRdeep_sRNAbench_intersect.bed \
  --sRNAbench-inter-table {RNA_MI_DIR}/sRNAbench_miRdeep_intersect.bed \
  --fc-mirdeep {SPECIES_DIR}/counts_sep/miRNA_mirdeep_counts.txt \
  --fc-pre-mirdeep {SPECIES_DIR}/counts_sep/miRNA_mirdeep_counts_flanked.txt \
  --fc-sRNAbench {SPECIES_DIR}/counts_sep/miRNA_sRNAbench_counts.txt \
  --fc-pre-sRNAbench {SPECIES_DIR}/counts_sep/miRNA_sRNAbench_counts_flanked.txt \
  -rm {SCRIPTS_DIR}/mirdeep_all_remaining_filtered.csv \
  -rs {SCRIPTS_DIR}/sRNAbench_all_remaining_filtered.csv \
  -l {LIBRARIES} \
  --sum-fc-thres 100
```

**Nematodes (non-Elegans)** — add BLAST:

```bash
  --blast-mirdeep {RNACENTRAL}/queries/{SPECIES}/miRdeep_blastn_compact \
  --blast-sRNAbench {RNACENTRAL}/queries/{SPECIES}/sRNAbench_blastn_compact
```

**Elegans only** — add miRBase/miRGeneDB BEDs, BLAST, `-mgff`, `--fc_mirbase`:

```bash
  --mirdeep-mibrase-inter {RNA_MI_DIR}/miRdeep_miRBase_intersect.bed \
  --mirdeep-mirgenedb-inter {RNA_MI_DIR}/miRdeep_miRGeneDB_intersect.bed \
  --sRNAbench-mibrase-inter {RNA_MI_DIR}/sRNAbench_miRBase_intersect.bed \
  --sRNAbench-mirgenedb-inter {RNA_MI_DIR}/sRNAbench_miRGeneDB_intersect.bed \
  --mirbase-mirgenedb-inter {RNA_MI_DIR}/miRBase_miRGeneDB_intersect.bed \
  --mirbase-mirdeep-inter {RNA_MI_DIR}/miRBase_miRdeep_intersect.bed \
  --mirbase-sRNAbench-inter {RNA_MI_DIR}/miRBase_sRNAbench_intersect.bed \
  --blast-mirdeep {RNACENTRAL}/queries/Elegans/miRdeep_blastn_compact \
  --blast-sRNAbench {RNACENTRAL}/queries/Elegans/sRNAbench_blastn_compact \
  --fc_mirbase {SPECIES_DIR}/counts_sep/miRNA_mirbase_counts.txt \
  -mgff {BASE}/mirbase_data/cel_mirbase_seq.gff3
```

Output: `intersections_table_{SPECIES}.xlsx` in `{RNA_MI_DIR}/` (filename matches `-s` exactly).

---

### 11. Downstream analysis

```bash
# All candidate FASTAs
python {REPO}/allCandidatesFasta.py \
  --all {RNA_MI_DIR}/intersections_table_{SPECIES}.xlsx \
  -s {SPECIES} {VARIANT}

# Structural features (Ziv)
python {REPO}/Ziv_feature_SOS.py \
  --precursors {RNA_MI_DIR}/all_candidates_hairpin.fasta \
  --mature {RNA_MI_DIR}/all_candidates_mature.fasta \
  --star {RNA_MI_DIR}/all_candidates_star.fasta \
  --species {SPECIES} {VARIANT} \
  --all-remaining {RNA_MI_DIR}/intersections_table_{SPECIES}.xlsx

# Statistics (+ Hofstenia genomic clustering)
python {REPO}/statistics.py \
  --all {BASE}/Ziv_Features/all_remaining_after_ziv_{SPECIES}.xlsx \
  -s {SPECIES} {VARIANT}
```

---

## Species-specific forks

| Fork | Species | Action |
|------|---------|--------|
| No BLAST | Hofstenia | Skip §9; omit `--blast-*` in §10 |
| miRBase + miRGeneDB | Elegans | Run `mirbaseToGFF3.py`; full intersection matrix; extra intersectionsTable args |
| Coordinate QC | Nematodes | Run `compare_genome_to_fasta.py` after unite pass 2; Elegans uses `{GENOME_FA_NO_WS}` |
| Hofstenia unite | Hofstenia | Use `--base-path {BASE}` on unite + processGoodCandidates; no `-seed` |
| Hofstenia reads | Hofstenia | STAR/sRNAbench input from `Fastq/.../filtered/`, not `TrimmedFastq/` |
| Elegans paths | Elegans | `Bash/` and `Genome/` (capital letters); others use lowercase `bash/`, `genome/` |
| `-s` casing | All | Always canonical `{SPECIES}`; lowercase fails `get_species_config()` |
| Ziv sheet | Nematodes vs Hofstenia | `(D) Structural Features` vs `(A) Unfiltered` |
| Sensitivity plots | Elegans | `run_sensitivity_plots: True` in config |

---

## Ziv structural thresholds

Build reference distributions once:

```bash
python {REPO}/mirgenedbThresholds.py
python {REPO}/Ziv_feature_SOS.py \
  --precursors {BASE}/mirgenedb_data_v3/ALL_mirgenedb_hairpin.fasta \
  --mature {BASE}/mirgenedb_data_v3/ALL_mirgenedb_mature.fasta \
  --star {BASE}/mirgenedb_data_v3/ALL_mirgenedb_star.fasta \
  --species miRGeneDB
python {REPO}/plot_series.py
```

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

## Optional steps

**mirTrace QC** (nematodes):

```bash
java -jar -Xms4G -Xmx4G {CONDA}/mirtrace.jar qc \
  --species cel --adapter AACTGTAGGCACCATCAAT --config config.txt
```

**Expression dynamics:**

```bash
python {REPO}/expression_dynamics.py \
  --all {BASE}/Ziv_Features/all_remaining_after_ziv_{SPECIES}.xlsx \
  --libraries {LIBRARIES} --time {TIMEPOINTS} -s {SPECIES}

python {REPO}/expression_dynamics_all.py \
  --all {BASE}/All_species/all_species_candidates.xlsx
```

**Cross-species seeds:**

```bash
cd {BASE}/All_species
python {REPO}/seed_frequency.py
```

**miRge** (after Ziv): regenerate FASTAs from the appropriate Ziv sheet → `create_combined_mature_star.py` → `generate_miRNA_GFF.py` → `miRge-build` → `mirge.sbatch` → `mirge_processing.py`. See species pipeline docs and `run_miRge.sh`.

---

## Alternate genome assemblies (`--variant new_genome`)

Reuse trimmed reads from the original track; rebuild indices on new scaffolds. Outputs under `{Species}_newGenome/`.

| Species | Genome FASTA |
|---------|--------------|
| Elegans | `Elegans_newGenome/genome/CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna` |
| Macrosperma | `Macrosperma_newGenome/genome/CMACR.caenorhabditis_macrosperma_JU2083_v2.scaffolds.fna` |
| Sulstoni | `Sulstoni_newGenome/genome/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna` |
| Hofstenia | `Hofstenia_newGenome/sRNA_PBonly/hofPB_v6.FINAL.fa` |

Add `{VARIANT}` = `--variant new_genome` to all config-aware Python scripts (or `-s {Species}_newGenome`).

---

## Library reference

### *C. elegans* (CE57–CE81)

Genome: `{GENOME_DIR}/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa`  
Whitespace-stripped (miRDeep): `{GENOME_DIR}/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa`  
Paths: `{BASH_DIR}` = `Elegans/Bash/`, `{GENOME_DIR}` = `Elegans/Genome/`, `{INDEX_BASENAME}` = `elegans`

### *C. macrosperma* (MR4–MR8)

Genome: `CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna`

### *C. sulstoni* (SR0–SR7)

Genome: `CSULS.caenorhabditis_sulstoni_JU2788_v1.scaffolds.fna`

### *Hofstenia miamia* (33 libraries)

AMP1–3, DI1–3, EC1–3, GA1–3, HL1–3, IST1–3, PDi1–3, PDii1–3, PH1–3, PL1–3, SMA1–3.  
Genome: `Hofstenia/Genome/refs/Hmia_ref/Hmia.030120.fasta`  
Reads: `Hofstenia/Fastq/Hmia_annotation/filtered/{LIBRARY}.filtered.fastq`  
`{INDEX_BASENAME}` = `hofstenia`

Full SRR ↔ library tables: `Pipeline Elegans.md`, `pipeline_config.py`.

---

## Adding a new species

1. Add entry to `SPECIES_CONFIG` in `pipeline_config.py`.
2. Place data under `{BASE}/{Species}/` following the directory layout.
3. Substitute template variables; loop over `cfg["libraries"]` for per-library steps.
4. Run phases 1–11; use Elegans as validation template if curated miRNAs exist in miRBase/miRGeneDB.

---

## Legacy doc drift (species-specific files)

Some commands in `Pipeline <Species>.md` predate strict `pipeline_config` validation and use **lowercase** `-s` (e.g. `-s elegans`). Current scripts require **canonical** `-s` (`Elegans`, `Macrosperma`, `Sulstoni`, `Hofstenia`). When assembling from templates, follow this file and `pipeline_config.py`, not legacy casing in species docs. Full expanded commands with correct paths are in [`Pipeline.md`](Pipeline.md).

---

## Citations

- cutadapt: DOI 10.14806/ej.17.1.200
- miRDeep2: Friedländer et al., NAR 2012 — https://doi.org/10.1093/nar/gkr688
- sRNAbench: Aparicio-Puerta et al., NAR 2019 — https://doi.org/10.1093/nar/gkz415
- bedtools: Quinlan & Hall, Bioinformatics 2010 — https://doi.org/10.1093/bioinformatics/btq033
- Nematode sRNA-seq: Nelson & Ambros, G3 2021 — PRJNA678899
