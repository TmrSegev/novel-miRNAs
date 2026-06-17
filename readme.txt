Order of running the files:

Elegans, Macrosperma and Sulstoni (per-library discovery; united downstream in scripts/):

Per library (in each sRNAbench output folder and mirdeep_out/<library>/):
nematodesRNAbenchFilter.py   (--filter-mc 100)
nematodeMirdeepFilter.py     (--filter-s 10 --exclude-c 1000 --filter-mc 100)

In <Species>/scripts/ (two-pass good_candidates, like Hofstenia):
nematodesRNAbenchGFF.py      (--goodcandidates False, then True)
process_debugging.py         (--tool sRNAbench)
nematodeMirdeepGFF.py        (--goodcandidates False, then True)
process_debugging.py         (--tool miRDeep)
compare_genome_to_fasta.py   (--mode discovery, after GFF pass 2)

overlapSenseAnti.py
mirbaseToGFF3.py (Elegans only)
filterSpacesBlastDB.py (needed only once, before creating the first blast DB)
add_flank_to_GFF.py
intersectionsTable.py        (with --fc-pre-mirdeep / --fc-pre-sRNAbench)
allCandidatesFasta.py
Ziv_feature_SOS.py (this script runs the Ziv_Git.py script)
statistics.py
expression_dynamics.py (Elegans only)

After analyzing all species:
seed_frequency.py
expression_dynamics_all.py

Legacy (single combined discovery run; superseded by scripts above):
sRNAbenchResultsToGFF3.py
mirdeepResultsToGFF3.py

________________________________
Hofstenia:
hofsteniaMirdeepFilter.py (run using sbatch file)
hofsteniaMirdeepGFF.py
hofsteniasRNAbenchFilter.py (run using sbatch file)
hofsteniasRNAbenchGFF.py
process_debugging_Hofstenia.py
overlapSenseAnti.py
add_flank_to_GFF.py
intersectionsTableHofstenia.py
allCandidatesFasta.py
Ziv_feature_SOS.py (this script runs the Ziv_Git.py script)
To calculate structural features thresholds:
	run mirgenedbThresholds.py
	then run Ziv_feature_SOS.py with species=miRGeneDB,
	and then run plot_series.py
statistics.py
