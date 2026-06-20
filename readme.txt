Order of running the files (unified pipeline — all species):

Per library (in each sRNAbench output folder and mirdeep_out/<library>/):
  srnabenchPerLibraryFilter.py   (default --filter-mc 10)
  mirdeepPerLibraryFilter.py     (default --filter-s 10 --exclude-c 100 --filter-mc 10)

In <Species>/scripts/ (two-pass good_candidates):
  srnabenchUniteGFF.py           (--goodcandidates False, then True; -seed optional, defaults from species config)
  processGoodCandidates.py       (--tool sRNAbench)
  mirdeepUniteGFF.py             (--goodcandidates False, then True)
  processGoodCandidates.py       (--tool miRDeep)
  compare_genome_to_fasta.py     (--mode discovery, after GFF pass 2; nematodes)

  overlapSenseAnti.py
  mirbaseToGFF3.py               (Elegans only)
  filterSpacesBlastDB.py         (nematodes with BLAST; once per DB)
  add_flank_to_GFF.py            (-s <Species> [--variant new_genome] [--base-path ...])
  intersectionsTable.py          (BLAST optional: required for nematodes, skipped for Hofstenia)
  allCandidatesFasta.py
  Ziv_feature_SOS.py
  statistics.py
  expression_dynamics.py         (Elegans only)

After analyzing all species:
  seed_frequency.py
  expression_dynamics_all.py

Seed file defaults (override with -seed):
  Hofstenia: mirbase_data/ALL_seed_family_from_mirgendb.csv (general miRs of all species)
  Nematodes (Elegans, Macrosperma, Sulstoni): mirbase_data/Seeds.txt (miRBase-based nematode miRNAs)

Alternate genome assembly (any species):
  Use --variant new_genome on config-aware scripts, or pass -s <Species>_newGenome as an alias.
  Directory layout: Charles_seq/<Species>_newGenome/ (mirdeep_out, scripts, STAR, …)
  GFF/fasta prefixes stay <Species>_* (same species label as the original assembly).
  Sequencing reads are shared with the original assembly track.

  Adding a new assembly later (e.g. Elegans, Sulstoni):
    1. Place scaffolds under Charles_seq/<Species>_newGenome/
    2. Add genome_fasta_subpath to NEW_GENOME_OVERRIDES in pipeline_config.py if not under genome/*.scaffolds.fna
    3. Copy/adapt bash from the existing species track; build bowtie/STAR/sRNAbench indices
    4. Run per-library discovery, then Python steps with --variant new_genome

Hofstenia notes:
  --variant new_genome selects the PacBio assembly (Hofstenia_newGenome/)

Macrosperma notes:
  --variant new_genome selects v2 scaffolds (Macrosperma_newGenome/); see docs/Pipeline Macrosperma.md

Backward-compatible wrappers (deprecated): nematode*GFF.py, hofstenia*GFF.py, nematode*Filter.py,
  hofstenia*Filter.py, process_debugging*.py, intersectionsTableHofstenia.py

Legacy (single combined discovery run; superseded):
  sRNAbenchResultsToGFF3.py
  mirdeepResultsToGFF3.py
