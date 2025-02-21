Order of running the files:

Elegans, Macrosperma and Sulstoni:
mirdeepResultsToGFF3.py
sRNAbenchResultsToGFF3.py
overlapSenseAnti.py
mirbaseToGFF3.py (Elegans only)
filterSpacesBlastDB.py (needed only once, before creating the first blast DB)
‏‏‏‏intersectionsTable.py
allCandidatesFasta.py
Ziv_feature_SOS.py (this script runs the Ziv_Git.py script)
statistics.py
expression_dynamics.py

After analyzing all species:
seed_frequency.py
expression_dynamics_all.py

________________________________
Hofstenia:
hofsteniaMirdeepFilter.py (run using sbatch file)
hofsteniaMirdeepGFF.py
hofsteniasRNAbenchFilter.py (run using sbatch file)
hofsteniasRNAbenchGFF.py
overlapSenseAnti.py
intersectionsTableHofstenia.py
allCandidatesFasta.py
Ziv_feature_SOS.py (this script runs the Ziv_Git.py script)
To calculate structural features thresholds:
	run mirgenedbThresholds.py
	then run Ziv_feature_SOS.py with species=miRGeneDB,
	and then run plot_series.py
statistics.py