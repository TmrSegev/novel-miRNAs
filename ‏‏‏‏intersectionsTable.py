import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

"""
Organize inputs
Create intersections table
Add BLAST results
Add feature counts
Add types
"""

# -----Getting Inputs-----
species = None
mirdeep_intersections_table_path = None
sRNAbench_intersections_table_path = None
blast_mirdeep_path = None
blast_sRNAbench_path = None
featurecounts_mirdeep_path = None
featurecounts_sRNAbench_path = None
remaining1_mirdeep_path = None
remaining2_mirdeep_path = None
remaining_sRNAbench_path = None
libraries = None
for i in range(1, len(sys.argv), 2):
    arg = sys.argv[i]
    if arg == '-s':
        species = sys.argv[i + 1]
    if arg == '--mirdeep-inter-table':
        mirdeep_intersections_table_path = sys.argv[i + 1]
    elif arg == '--sRNAbench-inter-table':
        sRNAbench_intersections_table_path = sys.argv[i + 1]
    elif arg == '--blast-mirdeep':
        blast_mirdeep_path = sys.argv[i + 1]
    elif arg == '--blast-sRNAbench':
        blast_sRNAbench_path = sys.argv[i + 1]
    elif arg == '--fc-mirdeep':
        featurecounts_mirdeep_path = sys.argv[i + 1]
    elif arg == '--fc-sRNAbench':
        featurecounts_sRNAbench_path = sys.argv[i + 1]
    elif arg == '-r1m':
        remaining1_mirdeep_path = sys.argv[i + 1]
    elif arg == '-r2m':
        remaining2_mirdeep_path = sys.argv[i + 1]
    elif arg == '-rs':
        remaining_sRNAbench_path = sys.argv[i + 1]
    elif arg == '-l':
        libraries = sys.argv[i + 1].split(',')
    elif arg == '--help' or arg == '-h':
        print(f'Manual:\n'
              f' -s <name>: name of species.\n'
              f' --mirdeep-inter-table <path>: path to bedtools -a mirdeep and -b sRNAbench intersection .bed file.\n'
              f' --sRNAbench-inter-table <path>: path to bedtools -a sRNAbench and -b mirdeep intersection .bed file.\n'
              f' --blast-mirdeep <path>: path to mirdeep blast results file.\n'
              f' --blast-sRNAbench <path>: path to sRNAbench blast results file.\n'
              f' --fc-mirdeep <path>: path to mirdeep featurecounts results file (full counts, not the summary file).\n'
              f' --fc-sRNAbench <path>: path to sRNAbench featurecounts results file (full counts, not the summary file).\n'
              f' -r1m <path>: path to the first remaining mirdeep candidates file, remaining_file_1.csv.\n'
              f' -r2m <path>: path to the second remaining mirdeep candidates file, remaining_file_2.csv.\n'
              f' -rs <path>: path to remaining sRNAbench candidates file, sRNAbench_remaining.csv.\n'
              f' -l <path>: list of sequencing libraries. Write the list seperated with commas, witout spaces. Example: library1,library2,library3 \n'
              )
        sys.exit()

# -----mirdeep intersections table:-----

mirdeep_intersections_table = pd.read_csv(mirdeep_intersections_table_path, sep='\t', names=['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2', 'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_sRNAbench', '.4', 'pre_miRNA2', 'Start_sRNAbench', 'End_sRNAbench', '.5', 'Strand_sRNAbench', '.6', 'Description_sRNAbench'])
mirdeep_intersections_table = mirdeep_intersections_table.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
mirdeep_intersections_table.to_csv("mirdeep_intersections_table", sep='\t')

# -----sRNAbench intersections table:-----

sRNAbench_intersections_table = pd.read_csv(sRNAbench_intersections_table_path, sep='\t', names=['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench', '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirdeep', '.4', 'pre_miRNA2', 'Start_mirdeep', 'End_mirdeep', '.5', 'Strand_mirdeep', '.6', 'Description_mirdeep'])
sRNAbench_intersections_table = sRNAbench_intersections_table.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
sRNAbench_intersections_table.to_csv("sRNAbench_intersections_table", sep='\t')

# -----Add blast results:-----
# ---miRdeep:
blast_mirdeep_orig = pd.read_csv(blast_mirdeep_path, sep='\t', names=['query_accession', 'subject_accession', '%_identical_matches', 'alignment_length', 'mismatches', 'gap_openings', 'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bitscore'])
blast_mirdeep_orig = blast_mirdeep_orig.drop_duplicates(subset=["query_accession"])

# Add strand column
mask = blast_mirdeep_orig['subject_start'] < blast_mirdeep_orig['subject_end']
blast_mirdeep_orig.loc[mask, 'strand'] = '+'
blast_mirdeep_orig['strand'].fillna('-', inplace=True)
blast_mirdeep = blast_mirdeep_orig.drop(['%_identical_matches', 'mismatches', 'gap_openings', 'subject_start', 'subject_end', 'bitscore'], axis=1)


# Create index column for blast
blast_mirdeep['index'] = blast_mirdeep['query_accession'].str.split('|')
blast_mirdeep['index'] = blast_mirdeep['index'].apply(lambda x : x[4])

# Create index column for mirdeep results
mirdeep_intersections_table['index'] = mirdeep_intersections_table['Description_mirdeep'].str.split(';')
mirdeep_intersections_table['index'] = mirdeep_intersections_table['index'].apply(lambda x : x[3])

# Merge mirdeep results and blast results
mirdeep_blast_intersections_table = pd.merge(mirdeep_intersections_table, blast_mirdeep, on='index', how='left')

# ---sRNAbench:

blast_sRNAbench_orig = pd.read_csv(blast_sRNAbench_path, sep='\t', names=['query_accession', 'subject_accession', '%_identical_matches', 'alignment_length', 'mismatches', 'gap_openings', 'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bitscore'])
blast_sRNAbench_orig = blast_sRNAbench_orig.drop_duplicates(subset=["query_accession"])

# Add strand column
mask = blast_sRNAbench_orig['subject_start'] < blast_sRNAbench_orig['subject_end']
blast_sRNAbench_orig.loc[mask, 'strand'] = '+'
blast_sRNAbench_orig['strand'].fillna('-', inplace=True)
blast_sRNAbench = blast_sRNAbench_orig.drop(['%_identical_matches', 'mismatches', 'gap_openings', 'subject_start', 'subject_end', 'bitscore'], axis=1)

# Create index column for blast
blast_sRNAbench['index'] = blast_sRNAbench['query_accession'].str.split('|')
blast_sRNAbench['index'] = blast_sRNAbench['index'].apply(lambda x : x[3])

# Create index column for sRNAbench results
sRNAbench_intersections_table['index'] = sRNAbench_intersections_table['Description_sRNAbench'].str.split(';')
sRNAbench_intersections_table['index'] = sRNAbench_intersections_table['index'].apply(lambda x : x[3])

# Merge sRNAbench results and blast results
sRNAbench_blast_intersections_table = pd.merge(sRNAbench_intersections_table, blast_sRNAbench, on='index', how='left')

# -----Add Featurecounts Results:-----
# ---miRDeep:
featurecounts_mirdeep = pd.read_csv(featurecounts_mirdeep_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
featurecounts_mirdeep = featurecounts_mirdeep.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
featurecounts_mirdeep = featurecounts_mirdeep.iloc[2:] # Drop the first 2 rows, which is readme info from featurecounts and not data


# Create mature/star and index column
featurecounts_mirdeep['index'] = featurecounts_mirdeep['Geneid'].str.split('|')
featurecounts_mirdeep['index'] = featurecounts_mirdeep['index'].apply(lambda x : x[4])

featurecounts_mirdeep['mature/star'] = featurecounts_mirdeep['Geneid'].str.split('|')
featurecounts_mirdeep['mature/star'] = featurecounts_mirdeep['mature/star'].apply(lambda x : x[2])

featurecounts_mirdeep = featurecounts_mirdeep.drop('Geneid', axis=1)

# Casting libraries columns to int64
cast_dict = {k: 'int64' for k in libraries}
featurecounts_mirdeep = featurecounts_mirdeep.astype(cast_dict)

# Separate df into mature and star
mature_counts = featurecounts_mirdeep[featurecounts_mirdeep['mature/star'] == 'm']
libraries_mature = [library + '_m' for library in libraries]
rename_dict = dict(zip(libraries, libraries_mature))
mature_counts = mature_counts.rename(columns=rename_dict)
mature_counts['sum_FC_m'] = np.zeros(len(mature_counts))
for library in libraries_mature:
    mature_counts['sum_FC_m'] += mature_counts[library]
mature_counts['sum_FC_m > 100?'] = np.where(mature_counts['sum_FC_m'] > 100, 1, 0)
mature_counts = mature_counts.drop('mature/star', axis=1)


star_counts = featurecounts_mirdeep[featurecounts_mirdeep['mature/star'] == 's']
libraries_star = [library + '_s' for library in libraries]
rename_dict = dict(zip(libraries, libraries_star))
star_counts = star_counts.rename(columns=rename_dict)
star_counts['sum_FC_s'] = np.zeros(len(star_counts))
for library in libraries_star:
    star_counts['sum_FC_s'] += star_counts[library]
star_counts['sum_FC_s > 100?'] = np.where(star_counts['sum_FC_s'] > 100, 1, 0)
star_counts = star_counts.drop('mature/star', axis=1)

# Merge mirdeep & blast results and featurecounts results
mirdeep_blast_m_intersections_table = pd.merge(mirdeep_blast_intersections_table, mature_counts, on='index', how='left')
mirdeep_blast_fc_intersections_table = pd.merge(mirdeep_blast_m_intersections_table, star_counts, on='index', how='left')
mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table.drop('index', axis=1)

# Extract readcounts columns
mirdeep_blast_fc_intersections_table['RC_m mirdeep'] = mirdeep_blast_fc_intersections_table["Description_mirdeep"].str.split(';', expand=True)[1]
mirdeep_blast_fc_intersections_table['RC_s mirdeep'] = mirdeep_blast_fc_intersections_table["Description_mirdeep"].str.split(';', expand=True)[2]
mirdeep_blast_fc_intersections_table['RC_m sRNAbench'] = mirdeep_blast_fc_intersections_table["Description_sRNAbench"].str.split(';', expand=True)[1]
mirdeep_blast_fc_intersections_table['RC_s sRNAbench'] = mirdeep_blast_fc_intersections_table["Description_sRNAbench"].str.split(';', expand=True)[2]
mirdeep_blast_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']] = mirdeep_blast_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']].fillna('0')
mirdeep_blast_fc_intersections_table['RC_m mirdeep'] = mirdeep_blast_fc_intersections_table['RC_m mirdeep'].str.replace('RC_m=', '').astype('int64')
mirdeep_blast_fc_intersections_table['RC_s mirdeep'] = mirdeep_blast_fc_intersections_table['RC_s mirdeep'].str.replace('RC_s=', '').astype('int64')
mirdeep_blast_fc_intersections_table['RC_m sRNAbench'] = mirdeep_blast_fc_intersections_table['RC_m sRNAbench'].str.replace('RC_m=', '').astype('int64')
mirdeep_blast_fc_intersections_table['RC_s sRNAbench'] = mirdeep_blast_fc_intersections_table['RC_s sRNAbench'].str.replace('RC_s=', '').astype('int64')

# Create diff columns
mirdeep_blast_fc_intersections_table['Diff Sum_FC_m / RC_m mirdeep'] = mirdeep_blast_fc_intersections_table['sum_FC_m'] / mirdeep_blast_fc_intersections_table['RC_m mirdeep']
mirdeep_blast_fc_intersections_table['Diff Sum_FC_m / RC_m sRNAbench'] = mirdeep_blast_fc_intersections_table['sum_FC_m'] / mirdeep_blast_fc_intersections_table['RC_m sRNAbench']
mirdeep_blast_fc_intersections_table['Diff Sum_FC_s / RC_s mirdeep'] = mirdeep_blast_fc_intersections_table['sum_FC_s'] / mirdeep_blast_fc_intersections_table['RC_s mirdeep']
mirdeep_blast_fc_intersections_table['Diff Sum_FC_s / RC_s sRNAbench'] = mirdeep_blast_fc_intersections_table['sum_FC_s'] / mirdeep_blast_fc_intersections_table['RC_s sRNAbench']

# Normalize featurecounts to reads per million
columns = [libraries_mature + libraries_star]
for column in columns:
    total = mirdeep_blast_fc_intersections_table[column].sum()
    mirdeep_blast_fc_intersections_table[column] = round((mirdeep_blast_fc_intersections_table[column] / total) * 1000000, 0)

# ADD SUM_FC_RPM!!!!

# ---sRNAbench:
featurecounts_sRNAbench = pd.read_csv(featurecounts_sRNAbench_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
featurecounts_sRNAbench = featurecounts_sRNAbench.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
featurecounts_sRNAbench = featurecounts_sRNAbench.iloc[2:] # Drop the first 2 rows, which is readme info from featurecounts and not data

# Create 5p/3p column
featurecounts_sRNAbench['5p/3p'] = featurecounts_sRNAbench['Geneid'].str.split('|', expand=True)[0]
featurecounts_sRNAbench['5p/3p'] = featurecounts_sRNAbench['5p/3p'].str.split('-')
featurecounts_sRNAbench['5p/3p'] = featurecounts_sRNAbench['5p/3p'].apply(lambda x : x[-1])
featurecounts_sRNAbench['5p/3p'] = featurecounts_sRNAbench['5p/3p'].str.split('_', expand=True)[0]
featurecounts_sRNAbench = featurecounts_sRNAbench.rename(columns={'5p/3p':'mature'})

# Create mature/star and index column
featurecounts_sRNAbench['index'] = featurecounts_sRNAbench['Geneid'].str.split('|')
featurecounts_sRNAbench['index'] = featurecounts_sRNAbench['index'].apply(lambda x : x[3])

featurecounts_sRNAbench['mature/star'] = featurecounts_sRNAbench['Geneid'].str.split('|')
featurecounts_sRNAbench['mature/star'] = featurecounts_sRNAbench['mature/star'].apply(lambda x : x[1])

featurecounts_sRNAbench = featurecounts_sRNAbench.drop('Geneid', axis=1)

# Casting libraries columns to int64
featurecounts_sRNAbench = featurecounts_sRNAbench.astype(cast_dict)

# Separate df into mature and star
mature_counts = featurecounts_sRNAbench[featurecounts_sRNAbench['mature/star'] == 'm']
rename_dict = dict(zip(libraries, libraries_mature))
mature_counts = mature_counts.rename(columns=rename_dict)
mature_counts['sum_FC_m'] = np.zeros(len(mature_counts))
for library in libraries_mature:
    mature_counts['sum_FC_m'] += mature_counts[library]
mature_counts['sum_FC_m > 100?'] = np.where(mature_counts['sum_FC_m'] > 100, 1, 0)
mature_counts = mature_counts.drop('mature/star', axis=1)

star_counts = featurecounts_sRNAbench[featurecounts_sRNAbench['mature/star'] == 's']
rename_dict = dict(zip(libraries, libraries_star))
star_counts = star_counts.rename(columns=rename_dict)
star_counts['sum_FC_s'] = np.zeros(len(star_counts))
for library in libraries_star:
    star_counts['sum_FC_s'] += star_counts[library]
star_counts['sum_FC_s > 100?'] = np.where(star_counts['sum_FC_s'] > 100, 1, 0)
star_counts = star_counts.drop(['mature/star', 'mature'], axis=1)

# Merge sRNAbench & blast results and featurecounts results
sRNAbench_blast_m_intersections_table = pd.merge(sRNAbench_blast_intersections_table, mature_counts, on='index', how='left')
sRNAbench_blast_fc_intersections_table = pd.merge(sRNAbench_blast_m_intersections_table, star_counts, on='index', how='left')
sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table.drop('index', axis=1)

# Extract readcounts columns
sRNAbench_blast_fc_intersections_table['RC_m mirdeep'] = sRNAbench_blast_fc_intersections_table["Description_mirdeep"].str.split(';', expand=True)[1]
sRNAbench_blast_fc_intersections_table['RC_s mirdeep'] = sRNAbench_blast_fc_intersections_table["Description_mirdeep"].str.split(';', expand=True)[2]
sRNAbench_blast_fc_intersections_table['RC_m sRNAbench'] = sRNAbench_blast_fc_intersections_table["Description_sRNAbench"].str.split(';', expand=True)[1]
sRNAbench_blast_fc_intersections_table['RC_s sRNAbench'] = sRNAbench_blast_fc_intersections_table["Description_sRNAbench"].str.split(';', expand=True)[2]
sRNAbench_blast_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']] = sRNAbench_blast_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']].fillna('0')
sRNAbench_blast_fc_intersections_table['RC_m mirdeep'] = sRNAbench_blast_fc_intersections_table['RC_m mirdeep'].str.replace('RC_m=', '').astype('int64')
sRNAbench_blast_fc_intersections_table['RC_s mirdeep'] = sRNAbench_blast_fc_intersections_table['RC_s mirdeep'].str.replace('RC_s=', '').astype('int64')
sRNAbench_blast_fc_intersections_table['RC_m sRNAbench'] = sRNAbench_blast_fc_intersections_table['RC_m sRNAbench'].str.replace('RC_m=', '').astype('int64')
sRNAbench_blast_fc_intersections_table['RC_s sRNAbench'] = sRNAbench_blast_fc_intersections_table['RC_s sRNAbench'].str.replace('RC_s=', '').astype('int64')

# Create diff columns
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_m / RC_m mirdeep'] = sRNAbench_blast_fc_intersections_table['sum_FC_m'] / sRNAbench_blast_fc_intersections_table['RC_m mirdeep']
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_m / RC_m sRNAbench'] = sRNAbench_blast_fc_intersections_table['sum_FC_m'] / sRNAbench_blast_fc_intersections_table['RC_m sRNAbench']
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_s / RC_s mirdeep'] = sRNAbench_blast_fc_intersections_table['sum_FC_s'] / sRNAbench_blast_fc_intersections_table['RC_s mirdeep']
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_s / RC_s sRNAbench'] = sRNAbench_blast_fc_intersections_table['sum_FC_s'] / sRNAbench_blast_fc_intersections_table['RC_s sRNAbench']

# -----Add Sequences:-----
# ---miRdeep:
remaining1_mirdeep = pd.read_csv(remaining1_mirdeep_path, sep='\t')
remaining2_mirdeep = pd.read_csv(remaining2_mirdeep_path, sep='\t')
remaining_mirdeep = pd.concat([remaining1_mirdeep, remaining2_mirdeep], ignore_index=True)
mirdeep_blast_fc_intersections_table['consensus mature sequence'] = remaining_mirdeep['consensus mature sequence'].str.upper()
mirdeep_blast_fc_intersections_table['consensus star sequence'] = remaining_mirdeep['consensus star sequence'].str.upper()
mirdeep_blast_fc_intersections_table['consensus precursor sequence'] = remaining_mirdeep['consensus precursor sequence'].str.upper()

# Create 5p/3p columns for mirdeep
def find_mature_index(row):
    index = row["consensus precursor sequence"].find(row["consensus mature sequence"])
    if index == 0:
        return '5p'
    elif index > 0:
        return '3p'
    elif index == -1:
        return 'error'

mirdeep_blast_fc_intersections_table['mature'] = mirdeep_blast_fc_intersections_table.apply(lambda row : find_mature_index(row), axis=1)

# Extract loop size
def loop_size_mirdeep(row):
    if (row["consensus precursor sequence"].find(str(row['consensus mature sequence'])) == -1) or (row["consensus precursor sequence"].find(str(row['consensus star sequence'])) == -1):
        return -1
    if row['mature'] == '5p':
        index_end_5p = len((str(row['consensus mature sequence'])))
        index_start_3p = row["consensus precursor sequence"].find(str(row['consensus star sequence']))
        loop_size = index_start_3p - index_end_5p
    elif row['mature'] == '3p':
        index_end_5p = len((str(row['consensus star sequence'])))
        index_start_3p = row["consensus precursor sequence"].find(str(row['consensus mature sequence']))
        loop_size = index_start_3p - index_end_5p
    return loop_size

mirdeep_blast_fc_intersections_table['loop_size'] = mirdeep_blast_fc_intersections_table.apply(lambda row : loop_size_mirdeep(row), axis=1)
# mirdeep_blast_fc_intersections_table['mature_size'] = np.where(mirdeep_blast_fc_intersections_table['mature'] == '5p', len(mirdeep_blast_fc_intersections_table['5pseq']), len(mirdeep_blast_fc_intersections_table['3pseq']))
mirdeep_blast_fc_intersections_table['mature_length'] = mirdeep_blast_fc_intersections_table['consensus mature sequence'].str.len()
mirdeep_blast_fc_intersections_table['star_length'] = mirdeep_blast_fc_intersections_table['consensus star sequence'].str.len()

#---sRNAbench:
remaining_sRNAbench = pd.read_csv(remaining_sRNAbench_path, sep='\t')
sRNAbench_blast_fc_intersections_table['5pseq'] = remaining_sRNAbench['5pseq'].str.replace('T', 'U')
sRNAbench_blast_fc_intersections_table['3pseq'] = remaining_sRNAbench['3pseq'].str.replace('T', 'U')
sRNAbench_blast_fc_intersections_table['hairpinSeq'] = remaining_sRNAbench['hairpinSeq'].str.replace('T', 'U')

#Removing flank and adjusting start/end
sRNAbench_blast_fc_intersections_table['hairpinSeq'] = sRNAbench_blast_fc_intersections_table['hairpinSeq'].str[11:-11]
sRNAbench_blast_fc_intersections_table['Start_sRNAbench'] = sRNAbench_blast_fc_intersections_table['Start_sRNAbench'] + 11
sRNAbench_blast_fc_intersections_table['End_sRNAbench'] = sRNAbench_blast_fc_intersections_table['End_sRNAbench'] - 11

# Extract loop size
def loop_size_sRNAbench(row):
    if (row["hairpinSeq"].find(str(row['5pseq'])) == -1) or (row["hairpinSeq"].find(str(row['3pseq'])) == -1):
        return -1
    index_end_5p = len((str(row['5pseq'])))
    index_start_3p = row["hairpinSeq"].find(str(row['3pseq']))
    loop_size = index_start_3p - index_end_5p
    return loop_size

sRNAbench_blast_fc_intersections_table['loop_size'] = sRNAbench_blast_fc_intersections_table.apply(lambda row : loop_size_sRNAbench(row), axis=1)

# -----Reorder columns:-----
mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table[['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep',
                                                                             'query_accession','subject_accession', 'alignment_length', 'query_start', 'query_end', 'e_value',
                                                                             'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench'] +
                                                                             libraries_mature + ['sum_FC_m', 'sum_FC_m > 100?', 'RC_m mirdeep', 'RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench'] +
                                                                             libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s mirdeep', 'RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench',
                                                                             'consensus mature sequence', 'consensus star sequence', 'consensus precursor sequence', 'mature', 'loop_size']
                                                                            ]

sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table[['Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench',
                                                                             'query_accession','subject_accession', 'alignment_length', 'query_start', 'query_end', 'e_value',
                                                                             'Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep'] +
                                                                             libraries_mature + ['sum_FC_m', 'sum_FC_m > 100?', 'RC_m sRNAbench', 'RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep'] +
                                                                             libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s sRNAbench', 'RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep',
                                                                             '5pseq', '3pseq', 'hairpinSeq', 'mature', 'loop_size']
                                                                                ]

# ------Add Types------
# ---miRdeep:
mirdeep_blast_fc_intersections_table['Type'] = np.where(mirdeep_blast_fc_intersections_table['Description_sRNAbench'] != '.', 1, 2)

# ---sRNAbench:
sRNAbench_blast_fc_intersections_table['Type'] = np.where(sRNAbench_blast_fc_intersections_table['Description_mirdeep'] != '.', 1, 3)

# ------Create unified sheet:------
unified = mirdeep_blast_fc_intersections_table.copy()

# Defining 5p/3p sequences in mirdeep
unified['5pseq'] = np.where(unified['mature'] == '5p', unified['consensus mature sequence'], unified['consensus star sequence'])
unified['3pseq'] = np.where(unified['mature'] == '3p', unified['consensus mature sequence'], unified['consensus star sequence'])

unified = unified.drop(['Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench', 'consensus mature sequence', 'consensus star sequence'], axis=1)

only_sRNAbench = sRNAbench_blast_fc_intersections_table[sRNAbench_blast_fc_intersections_table['Type'] == 3]
only_sRNAbench = only_sRNAbench.drop(['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep'], axis=1)
unified = unified.rename(columns={'Chr_mirdeep':'Chr', 'Start_mirdeep':'Start', 'End_mirdeep':'End', 'Strand_mirdeep':'Strand', 'Description_mirdeep':'Description', 'consensus precursor sequence':'hairpinSeq'})
only_sRNAbench = only_sRNAbench.rename(columns={'Chr_sRNAbench':'Chr', 'Start_sRNAbench':'Start', 'End_sRNAbench':'End', 'Strand_sRNAbench':'Strand', 'Description_sRNAbench':'Description'})
unified = pd.concat([unified, only_sRNAbench], axis=0, ignore_index=True)
# reorder sequences columns
columns = list(unified.columns)
i1, i2, i3 = columns.index('mature'), columns.index('5pseq'), columns.index('3pseq')
columns[i1], columns[i2], columns[i3] = columns[i2], columns[i3], columns[i1]
unified = unified[columns]

# # ---Normalize featurecounts to reads per million
# columns = [libraries_mature + ['sum_FC_m'] + libraries_star + ['sum_FC_s']]
# for column in columns:
#     total = unified[column].sum()
#     unified[column] = round((unified[column] / total) * 1000000, 0)

# --- Extract seed
unified['Seed'] = unified['Description'].str.split(';', expand=True)[4]
unified['Seed'] = unified['Seed'].str.split('|', expand=True)[0] # Remove sense / antisense / overlap

seed_families = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/cel-mirgenedb-families.csv')
seed_families = seed_families[['MiRBase ID', 'Family', 'Seed']]
# seed_families.loc[len(seed_families)] = ['test-mir-71', 'TEST', 'UACCUUG']
unified = pd.merge(unified, seed_families, left_on='Seed', right_on='MiRBase ID', how='left')
unified = pd.merge(unified, seed_families, left_on='Seed_x', right_on='Seed', how='left')
# unified = unified.drop_duplicates(subset='Description')
unified['Family_x'].fillna(' ', inplace=True)
unified['Family_y'].fillna(' ', inplace=True)
unified['Family'] = unified['Family_x'] + unified['Family_y']
unified['Family'] = unified['Family'].str.replace('  ', 'UNKNOWN')

unified = unified.drop(['Seed', 'MiRBase ID_x', 'Family_x', 'Family_y'], axis=1)
unified = unified.rename(columns={'Seed_x':'Candidate seed/ID', 'Seed_y':'Seed_mirGeneDB', 'MiRBase ID_y':'MiRBase ID_mirGeneDB'})
unified = unified.dropna(axis=1, thresh=1) # drop empty columns if there are any
unified['Seed'] = np.where(unified["mature"] == '5p', unified["5pseq"].str[1:8], unified["3pseq"].str[1:8])
unified = unified.reindex(columns=[col for col in unified.columns if col != 'Type'] + ['Type']) # move 'type' column to last position

# --- Removing duplicates
mask = unified.duplicated(subset=unified.columns.drop('Seed_mirGeneDB'), keep=False)
unified = unified.loc[~mask | (unified['Seed'] == unified['Seed_mirGeneDB'])]

# ---Creating a families by type pivot table
families_by_type = pd.pivot_table(unified, values='Description', index='Family', columns='Type', aggfunc='count')
families_by_type.fillna(0, inplace=True)
families_by_type.drop("UNKNOWN", axis=0, inplace=True)
families_by_type.plot(kind='bar', figsize=(14, 10), stacked=True, title="{} counts of known families in each type".format(species))
plt.ylabel("Counts")
plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
# plt.show()
plt.savefig("{}_family_counts_per_type.png".format(species))

# # for every type, create a bar chart of all the families
# types = unified['Type'].unique()
# for type in types:
#     filtered_unified = unified[unified['Type'] == type]
#     print(filtered_unified['Family'].value_counts())
#     filtered_unified['Family'].value_counts().plot(kind='bar', figsize=(15, 10))
#     plt.xticks(rotation=45, ha='right')
#     plt.title("Counts by families for candidates of type " + str(type))
#     plt.show()

# -----Save each intersections table as a sheet in one excel file:-----

writer = pd.ExcelWriter('intersections_table_{}.xlsx'.format(species))
mirdeep_blast_fc_intersections_table.to_excel(writer, sheet_name='miRdeep')
sRNAbench_blast_fc_intersections_table.to_excel(writer, sheet_name='sRNAbench')
unified.to_excel(writer, sheet_name='all_candidates')
blast_mirdeep_orig.to_excel(writer, sheet_name='blast_miRdeep')
blast_sRNAbench_orig.to_excel(writer, sheet_name='blast_sRNAbench')
writer.save()

