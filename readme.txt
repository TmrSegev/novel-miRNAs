Order of running the files (unified pipeline — all species):

Per library (in each sRNAbench output folder and mirdeep_out/<library>/):
  srnabenchPerLibraryFilter.py   (default --filter-mc 10; nematodes historically: 100)
  mirdeepPerLibraryFilter.py     (default --filter-s 10 --exclude-c 100 --filter-mc 10;
                                  nematodes historically: --filter-mc 100 --exclude-c 1000)

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
  Nematodes (Elegans, Macrosperma, Sulstoni): mirbase_data/Seeds.txt
  Hofstenia: mirbase_data/ALL_seed_family_from_mirgendb.csv

Hofstenia notes:
  Use --variant new_genome for PacBio genome track

Backward-compatible wrappers (deprecated): nematode*GFF.py, hofstenia*GFF.py, nematode*Filter.py,
  hofstenia*Filter.py, process_debugging*.py, intersectionsTableHofstenia.py

Legacy (single combined discovery run; superseded):
  sRNAbenchResultsToGFF3.py
  mirdeepResultsToGFF3.py
