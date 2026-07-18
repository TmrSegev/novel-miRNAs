Scripts directory: /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/
Invoke all pipeline scripts by absolute path from that directory (do not copy scripts into species folders).
Pipeline documentation:
  docs/pipeline_v3.md          — canonical workflow (template reference, Phases 1–13)
  docs/Pipeline <Species>.md — species-specific notes and expanded examples
Legacy duplicate copies outside the repo were renamed to *_archived.py; see docs/archived_scripts_manifest.txt.

Order of running the files (see docs/pipeline_v3.md):

Phase 2 — per library:
  srnabenchPerLibraryFilter.py   (default --filter-mc 10)
  mirdeepPerLibraryFilter.py     (default --filter-s 10 --exclude-c 100 --filter-mc 10)

Phase 3 — curation in <Species>/scripts/ (good_candidates: steps A → B → C):
  srnabenchUniteGFF.py           (step A: --goodcandidates False; step C: True)
  processGoodCandidates.py       (step B; --tool sRNAbench or miRDeep)
  mirdeepUniteGFF.py             (step A / step C)
  compare_genome_to_fasta.py     (--mode discovery; nematodes)
  overlapSenseAnti.py

Phase 4 — quantification:
  add_flank_to_GFF.py            (-s <Species> [--variant new_genome])
  filterSpacesBlastDB.py         (nematodes; once per DB)
  STAR/featureCounts + blastn     (see docs/pipeline_v3.md Phases 7–8)

Phase 5 — integration:
  mirbaseToGFF3.py               (Elegans only)
  bedtools cross-intersections + intersectionsTable.py

Phase 6 — final filters:
  allCandidatesFasta.py
  Ziv_feature_SOS.py             (structural filter; input = intersections table)
  statistics.py

Optional: expression_dynamics.py (Elegans, Hofstenia)

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
  --variant new_genome selects v2 scaffolds (Macrosperma_newGenome/); see docs/pipeline_v3.md

Backward-compatible wrappers (deprecated): nematode*GFF.py, hofstenia*GFF.py, nematode*Filter.py,
  hofstenia*Filter.py, process_debugging*.py, intersectionsTableHofstenia.py

Legacy (single combined discovery run; superseded):
  sRNAbenchResultsToGFF3.py
  mirdeepResultsToGFF3.py
