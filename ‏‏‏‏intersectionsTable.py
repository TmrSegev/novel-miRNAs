import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

"""
GETTING INPUTS
CREATE INTERSECTIONS TABLES
ADD BLAST RESULTS
ADD FEATURE COUNTS
ADD SEQUENCES
REORDER COLUMNS
ADD TYPES
CREATE ALL CANDIDATES SHEET
STATISTICS
SAVE TO EXCEL
"""

# -----GETTING INPUTS-----
species = None
mirdeep_intersections_table_path = None
mirdeep_mibrase_inter = None
mirdeep_mirgenedb_inter = None
sRNAbench_intersections_table_path = None
sRNAbench_mibrase_inter = None
sRNAbench_mirgenedb_inter = None
mirbase_mirgenedb_inter = None
mirbase_mirdeep_inter = None
mirbase_sRNAbench_inter = None
blast_mirdeep_path = None
blast_sRNAbench_path = None
featurecounts_mirdeep_path = None
featurecounts_sRNAbench_path = None
featurecounts_mirbase_path = None
remaining1_mirdeep_path = None
remaining2_mirdeep_path = None
remaining_sRNAbench_path = None
libraries = None
sum_fc_thres = 100
for i in range(1, len(sys.argv), 2):
    arg = sys.argv[i]
    if arg == '-s':
        species = sys.argv[i + 1]
    elif arg == '--mirdeep-inter-table':
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
    elif arg == '--fc_mirbase':
        featurecounts_mirbase_path = sys.argv[i + 1]
    elif arg == '-r1m':
        remaining1_mirdeep_path = sys.argv[i + 1]
    elif arg == '-r2m':
        remaining2_mirdeep_path = sys.argv[i + 1]
    elif arg == '-rs':
        remaining_sRNAbench_path = sys.argv[i + 1]
    elif arg == '-l':
        libraries = sys.argv[i + 1].split(',')
    elif arg == '--sum-fc-thres':
        sum_fc_thres = int(sys.argv[i + 1])
    elif arg == '--mirdeep-mibrase-inter':
        mirdeep_mibrase_inter = sys.argv[i + 1]
    elif arg == '--mirdeep-mirgenedb-inter':
        mirdeep_mirgenedb_inter = sys.argv[i + 1]
    elif arg == '--sRNAbench-mibrase-inter':
        sRNAbench_mibrase_inter = sys.argv[i + 1]
    elif arg == '--sRNAbench-mirgenedb-inter':
        sRNAbench_mirgenedb_inter = sys.argv[i + 1]
    elif arg == '--mirbase-mirgenedb-inter':
        mirbase_mirgenedb_inter = sys.argv[i + 1]
    elif arg == '--mirbase-mirdeep-inter':
        mirbase_mirdeep_inter = sys.argv[i + 1]
    elif arg == '--mirbase-sRNAbench-inter':
        mirbase_sRNAbench_inter = sys.argv[i + 1]
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
              f' -l <list>: list of sequencing libraries. Write the list seperated with commas, witout spaces. Example: library1,library2,library3 \n'
              f' --sum-fc-thres <int>: filtering threshold. Any candidates with sum_fc_m <= threshold will be filtered.\n'
              f'\nElegans only parameters:\n'
              f' --mirdeep-mibrase-inter <path>: path to bedtools -a mirdeep and -b mirbase intersection .bed file.\n'
              f' --mirdeep-mirgenedb-inter <path>: path to bedtools -a mirdeep and -b mirgenedb intersection .bed file.\n'
              f' --sRNAbench-mibrase-inter <path>: path to bedtools -a sRNAbench and -b mirbase intersection .bed file.\n'
              f' --sRNAbench-mirgenedb-inter <path>: path to bedtools -a sRNAbench and -b mirgenedb intersection .bed file.\n'
              f' --mirbase-mirgenedb-inter <path>: path to bedtools -a mirbase and -b mirgenedb intersection .bed file.\n'
              f' --mirbase-mirdeep-inter <path>: path to bedtools -a mirbase and -b mirdeep intersection .bed file.\n'
              f' --mirbase-sRNAbench-inter <path>: path to bedtools -a mirbase and -b sRNAbench intersection .bed file.\n'
              f' --fc-mirbase <path>: path to mirbase featurecounts results file (full counts, not the summary file).\n'
              )
        sys.exit()

# -----CREATE INTERSECTIONS TABLES-----
# -----mirdeep intersections table:-----

mirdeep_intersections_table = pd.read_csv(mirdeep_intersections_table_path, sep='\t', names=['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2', 'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_sRNAbench', '.4', 'pre_miRNA2', 'Start_sRNAbench', 'End_sRNAbench', '.5', 'Strand_sRNAbench', '.6', 'Description_sRNAbench'])
mirdeep_intersections_table = mirdeep_intersections_table.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
# mirdeep_intersections_table.to_csv("mirdeep_intersections_table", sep='\t')

if (species == 'Elegans') or (species == 'elegans'):
    mirdeep_mirbase = pd.read_csv(mirdeep_mibrase_inter, sep='\t',
                                  names=['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2',
                                         'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_mirbase', '.4',
                                         'pre_miRNA2', 'Start_mirbase', 'End_mirbase', '.5', 'Strand_mirbase', '.6',
                                         'Description_mirbase'])
    mirdeep_mirgenedb = pd.read_csv(mirdeep_mirgenedb_inter, sep='\t',
                                    names=['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2',
                                           'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_mirgenedb', '.4',
                                           'pre_miRNA2', 'Start_mirgenedb', 'End_mirgenedb', '.5', 'Strand_mirgenedb',
                                           '.6', 'Description_mirgenedb'])

    mirdeep_mirbase = mirdeep_mirbase.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
    mirdeep_mirgenedb = mirdeep_mirgenedb.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)

    mirdeep_intersections_table['T/F_sRNAbench'] = (mirdeep_intersections_table['Description_sRNAbench'] != '.').astype(
        int)  # Used for classifying types

    mirdeep_sRNAbench_mirbase = pd.merge(mirdeep_intersections_table, mirdeep_mirbase.iloc[:, 4:10], on='Description_mirdeep',
                                         how='left')
    mirdeep_sRNAbench_mirbase['T/F_mirbase'] = (mirdeep_sRNAbench_mirbase['Description_mirbase'] != '.').astype(
        int)  # Used for classifying types
    mirdeep_intersections_table = pd.merge(mirdeep_sRNAbench_mirbase, mirdeep_mirgenedb.iloc[:, 4:10],
                                           on='Description_mirdeep', how='left')
    mirdeep_intersections_table['T/F_mirgenedb'] = (mirdeep_intersections_table['Description_mirgenedb'] != '.').astype(
        int)  # Used for classifying types
    # mirdeep_intersections_table.to_csv("mirdeep_intersections_table", sep='\t')


# -----sRNAbench intersections table:-----

sRNAbench_intersections_table = pd.read_csv(sRNAbench_intersections_table_path, sep='\t', names=['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench', '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirdeep', '.4', 'pre_miRNA2', 'Start_mirdeep', 'End_mirdeep', '.5', 'Strand_mirdeep', '.6', 'Description_mirdeep'])
sRNAbench_intersections_table = sRNAbench_intersections_table.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
# sRNAbench_intersections_table.to_csv("sRNAbench_intersections_table", sep='\t')

if (species == 'Elegans') or (species == 'elegans'):
    sRNAbench_mirbase = pd.read_csv(sRNAbench_mibrase_inter, sep='\t',
                                    names=['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench',
                                           '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirbase', '.4',
                                           'pre_miRNA2', 'Start_mirbase', 'End_mirbase', '.5', 'Strand_mirbase', '.6',
                                           'Description_mirbase'])
    sRNAbench_mirgenedb = pd.read_csv(sRNAbench_mirgenedb_inter, sep='\t',
                                      names=['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench',
                                             '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirgenedb',
                                             '.4', 'pre_miRNA2', 'Start_mirgenedb', 'End_mirgenedb', '.5',
                                             'Strand_mirgenedb', '.6', 'Description_mirgenedb'])

    sRNAbench_mirbase = sRNAbench_mirbase.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
    sRNAbench_mirgenedb = sRNAbench_mirgenedb.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'],
                                                   axis=1)

    sRNAbench_intersections_table['T/F_mirdeep'] = (sRNAbench_intersections_table['Description_mirdeep'] != '.').astype(
        int)  # Used for classifying types

    sRNAbench_mirdeep_mirbase = pd.merge(sRNAbench_intersections_table, sRNAbench_mirbase.iloc[:, 4:10], on='Description_sRNAbench',
                                         how='left')
    sRNAbench_mirdeep_mirbase['T/F_mirbase'] = (sRNAbench_mirdeep_mirbase['Description_mirbase'] != '.').astype(
        int)  # Used for classifying types
    sRNAbench_intersections_table = pd.merge(sRNAbench_mirdeep_mirbase, sRNAbench_mirgenedb.iloc[:, 4:10],
                                             on='Description_sRNAbench', how='left')
    sRNAbench_intersections_table['T/F_mirgenedb'] = (
                sRNAbench_intersections_table['Description_mirgenedb'] != '.').astype(int)  # Used for classifying types
    # sRNAbench_intersections_table.to_csv("sRNAbench_intersections_table", sep='\t')

if (species == 'Elegans') or (species == 'elegans'):
    # -----mirbase intersections table:-----
    mirbase_mirgenedb = pd.read_csv(mirbase_mirgenedb_inter, sep='\t',
                                    names=['Chr_mirbase', '.1', 'pre_miRNA1', 'Start_mirbase', 'End_mirbase', '.2',
                                           'Strand_mirbase', '.3', 'Description_mirbase', 'Chr_mirgenedb', '.4',
                                           'pre_miRNA2', 'Start_mirgenedb', 'End_mirgenedb', '.5', 'Strand_mirgenedb',
                                           '.6', 'Description_mirgenedb'])
    mirbase_mirdeep = pd.read_csv(mirbase_mirdeep_inter, sep='\t',
                                  names=['Chr_mirbase', '.1', 'pre_miRNA1', 'Start_mirbase', 'End_mirbase', '.2',
                                         'Strand_mirbase', '.3', 'Description_mirbase', 'Chr_mirdeep', '.4',
                                         'pre_miRNA2', 'Start_mirdeep', 'End_mirdeep', '.5', 'Strand_mirdeep', '.6',
                                         'Description_mirdeep'])
    mirbase_sRNAbench = pd.read_csv(mirbase_sRNAbench_inter, sep='\t',
                                    names=['Chr_mirbase', '.1', 'pre_miRNA1', 'Start_mirbase', 'End_mirbase', '.2',
                                           'Strand_mirbase', '.3', 'Description_mirbase', 'Chr_sRNAbench', '.4',
                                           'pre_miRNA2', 'Start_sRNAbench', 'End_sRNAbench', '.5', 'Strand_sRNAbench',
                                           '.6', 'Description_sRNAbench'])

    mirbase_mirgenedb = mirbase_mirgenedb.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
    mirbase_mirdeep = mirbase_mirdeep.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
    mirbase_sRNAbench = mirbase_sRNAbench.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)

    mirbase_mirgenedb['T/F_mirgenedb'] = (mirbase_mirgenedb['Description_mirgenedb'] != '.').astype(
        int)  # Used for classifying types

    mirbase_mirgenedb_mirdeep = pd.merge(mirbase_mirgenedb, mirbase_mirdeep.iloc[:, 4:10], on='Description_mirbase',
                                         how='left')
    mirbase_mirgenedb_mirdeep['T/F_mirdeep'] = (mirbase_mirgenedb_mirdeep['Description_mirdeep'] != '.').astype(
        int)  # Used for classifying types
    mirbase_intersections_table = pd.merge(mirbase_mirgenedb_mirdeep, mirbase_sRNAbench.iloc[:, 4:10],
                                           on='Description_mirbase', how='left')
    mirbase_intersections_table['T/F_sRNAbench'] = (mirbase_intersections_table['Description_sRNAbench'] != '.').astype(
        int)  # Used for classifying types
    # mirbase_intersections_table.to_csv("mirbase_intersections_table", sep='\t')

# -----ADD BLAST RESULTS-----
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

# -----ADD FEATURE COUNTS:-----
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
# mature_counts['sum_FC_m'] = np.zeros(len(mature_counts))
# for library in libraries_mature:
#     mature_counts['sum_FC_m'] += mature_counts[library]
mature_counts['sum_FC_m'] = mature_counts[libraries_mature].sum(axis=1)
mature_counts = mature_counts.drop('mature/star', axis=1)


star_counts = featurecounts_mirdeep[featurecounts_mirdeep['mature/star'] == 's']
libraries_star = [library + '_s' for library in libraries]
rename_dict = dict(zip(libraries, libraries_star))
star_counts = star_counts.rename(columns=rename_dict)
# star_counts['sum_FC_s'] = np.zeros(len(star_counts))
# for library in libraries_star:
#     star_counts['sum_FC_s'] += star_counts[library]
star_counts['sum_FC_s'] = star_counts[libraries_star].sum(axis=1)
star_counts['sum_FC_s > 100?'] = np.where(star_counts['sum_FC_s'] > 100, 1, 0)
star_counts = star_counts.drop('mature/star', axis=1)

# Merge mirdeep & blast results and featurecounts results
mirdeep_blast_m_intersections_table = pd.merge(mirdeep_blast_intersections_table, mature_counts, on='index', how='left')
mirdeep_blast_fc_intersections_table = pd.merge(mirdeep_blast_m_intersections_table, star_counts, on='index', how='left')
mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table.drop('index', axis=1)

# filter by sum_fc_m < threshold
mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table[mirdeep_blast_fc_intersections_table['sum_FC_m'] > sum_fc_thres]

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
total = mirdeep_blast_fc_intersections_table[libraries_mature + libraries_star].sum()
mature_rpm = [column + "_rpm" for column in libraries_mature]
star_rpm = [column + "_rpm" for column in libraries_star]
mirdeep_blast_fc_intersections_table[mature_rpm + star_rpm] = round((mirdeep_blast_fc_intersections_table[libraries_mature + libraries_star] / total) * 1000000, 0)

mirdeep_blast_fc_intersections_table['sum_FC_m_rpm'] = np.zeros(len(mirdeep_blast_fc_intersections_table))
for library in mature_rpm:
    mirdeep_blast_fc_intersections_table['sum_FC_m_rpm'] += mirdeep_blast_fc_intersections_table[library]
mirdeep_blast_fc_intersections_table['sum_FC_s_rpm'] = np.zeros(len(mirdeep_blast_fc_intersections_table))
for library in star_rpm:
    mirdeep_blast_fc_intersections_table['sum_FC_s_rpm'] += mirdeep_blast_fc_intersections_table[library]
mirdeep_blast_fc_intersections_table['mean_m_rpm'] = round(mirdeep_blast_fc_intersections_table[mature_rpm].mean(axis=1), 2)
mirdeep_blast_fc_intersections_table['mean_s_rpm'] = round(mirdeep_blast_fc_intersections_table[star_rpm].mean(axis=1), 2)

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
# mature_counts['sum_FC_m'] = np.zeros(len(mature_counts))
# for library in libraries_mature:
#     mature_counts['sum_FC_m'] += mature_counts[library]
mature_counts['sum_FC_m'] = mature_counts[libraries_mature].sum(axis=1)
mature_counts = mature_counts.drop('mature/star', axis=1)

star_counts = featurecounts_sRNAbench[featurecounts_sRNAbench['mature/star'] == 's']
rename_dict = dict(zip(libraries, libraries_star))
star_counts = star_counts.rename(columns=rename_dict)
# star_counts['sum_FC_s'] = np.zeros(len(star_counts))
# for library in libraries_star:
#     star_counts['sum_FC_s'] += star_counts[library]
star_counts['sum_FC_s'] = star_counts[libraries_star].sum(axis=1)
star_counts['sum_FC_s > 100?'] = np.where(star_counts['sum_FC_s'] > 100, 1, 0)
star_counts = star_counts.drop(['mature/star', 'mature'], axis=1)

# Merge sRNAbench & blast results and featurecounts results
sRNAbench_blast_m_intersections_table = pd.merge(sRNAbench_blast_intersections_table, mature_counts, on='index', how='left')
sRNAbench_blast_fc_intersections_table = pd.merge(sRNAbench_blast_m_intersections_table, star_counts, on='index', how='left')
sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table.drop('index', axis=1)

# filter by sum_fc_m < threshold
sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table[sRNAbench_blast_fc_intersections_table['sum_FC_m'] > sum_fc_thres]

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

# Normalize featurecounts to reads per million
total = sRNAbench_blast_fc_intersections_table[libraries_mature + libraries_star].sum()
sRNAbench_blast_fc_intersections_table[mature_rpm + star_rpm] = round((sRNAbench_blast_fc_intersections_table[libraries_mature + libraries_star] / total) * 1000000, 0)

sRNAbench_blast_fc_intersections_table['sum_FC_m_rpm'] = np.zeros(len(sRNAbench_blast_fc_intersections_table))
for library in mature_rpm:
    sRNAbench_blast_fc_intersections_table['sum_FC_m_rpm'] += sRNAbench_blast_fc_intersections_table[library]
sRNAbench_blast_fc_intersections_table['sum_FC_s_rpm'] = np.zeros(len(sRNAbench_blast_fc_intersections_table))
for library in star_rpm:
    sRNAbench_blast_fc_intersections_table['sum_FC_s_rpm'] += sRNAbench_blast_fc_intersections_table[library]
sRNAbench_blast_fc_intersections_table['mean_m_rpm'] = round(sRNAbench_blast_fc_intersections_table[mature_rpm].mean(axis=1), 2)
sRNAbench_blast_fc_intersections_table['mean_s_rpm'] = round(sRNAbench_blast_fc_intersections_table[star_rpm].mean(axis=1), 2)

if (species == 'Elegans') or (species == 'elegans'):
    # -----miRBase:
    featurecounts_mirbase = pd.read_csv(featurecounts_mirbase_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
    featurecounts_mirbase = featurecounts_mirbase.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
    featurecounts_mirbase = featurecounts_mirbase.iloc[2:] # Drop the first 2 rows, which is readme info from featurecounts and not data


    # Create index column for featurecounts
    featurecounts_mirbase['index'] = featurecounts_mirbase['Geneid'].str.split('|')
    featurecounts_mirbase['index'] = featurecounts_mirbase['index'].apply(lambda x : x[3])
    featurecounts_mirbase['index'] = featurecounts_mirbase['index'].str.replace('Derives_from=MI', '')

    # Create 5p/3p columns
    featurecounts_mirbase['5p/3p'] = featurecounts_mirbase['Geneid'].str.split('|')
    featurecounts_mirbase['5p/3p'] = featurecounts_mirbase['5p/3p'].apply(lambda x : x[2])
    featurecounts_mirbase['5p/3p'] = featurecounts_mirbase['5p/3p'].str[-2:]
    featurecounts_mirbase = featurecounts_mirbase.drop('Geneid', axis=1)

    # Casting libraries columns to int64
    featurecounts_mirbase = featurecounts_mirbase.astype(cast_dict)

    # Separate df into 5p and 3p
    counts_mb_5p = featurecounts_mirbase[featurecounts_mirbase['5p/3p'] == '5p']
    libraries_5p = [library + '_5p' for library in libraries]
    rename_dict_5p = dict(zip(libraries, libraries_5p))
    counts_mb_5p = counts_mb_5p.rename(columns=rename_dict_5p)
    counts_mb_5p = counts_mb_5p.drop('5p/3p', axis=1)

    counts_mb_3p = featurecounts_mirbase[featurecounts_mirbase['5p/3p'] == '3p']
    libraries_3p = [library + '_3p' for library in libraries]
    rename_dict_3p = dict(zip(libraries, libraries_3p))
    counts_mb_3p = counts_mb_3p.rename(columns=rename_dict_3p)
    counts_mb_3p = counts_mb_3p.drop('5p/3p', axis=1)

    counts_no_5p3p = featurecounts_mirbase[(featurecounts_mirbase['5p/3p'] != '3p') & (featurecounts_mirbase['5p/3p'] != '5p')]

    # sum 5p and 3p to determine mature/star
    numeric_counts_mb_5p = counts_mb_5p[libraries_5p].astype(int)
    counts_mb_5p['sum'] = numeric_counts_mb_5p.sum(axis=1)
    numeric_counts_mb_3p = counts_mb_3p[libraries_3p].astype(int)
    counts_mb_3p['sum'] = numeric_counts_mb_3p.sum(axis=1)

    print(counts_mb_5p['sum'])

    # Create empty mature df and star df, iterate the rows of 5p and 3p and add to the relavant.
    mature_df = pd.DataFrame(columns=libraries_mature)
    star_df = pd.DataFrame(columns=libraries_star)
    for i in range(0, len(counts_mb_5p)):
        row_5p = counts_mb_5p.iloc[i]
        row_3p = counts_mb_3p.iloc[i]
        if row_5p['sum'] > row_3p['sum']: # if 5p is the mature
            rename_dict = dict(zip(libraries_5p, libraries_mature))
            row_5p = row_5p.rename(rename_dict)
            rename_dict = dict(zip(libraries_3p, libraries_star))
            row_3p = row_3p.rename(rename_dict)
            mature_df = mature_df.append(row_5p)
            star_df = star_df.append(row_3p)
        else: # else 3p is the mature
            rename_dict = dict(zip(libraries_5p, libraries_star))
            row_5p = row_5p.rename(rename_dict)
            rename_dict = dict(zip(libraries_3p, libraries_mature))
            row_3p = row_3p.rename(rename_dict)
            mature_df = mature_df.append(row_3p)
            star_df = star_df.append(row_5p)

    mature_df = mature_df.rename(columns={'sum':'sum_FC_m'})
    star_df = star_df.rename(columns={'sum':'sum_FC_s'})
    star_df['sum_FC_s > 100?'] = np.where(star_df['sum_FC_s'] > 100, 1, 0)

    # Those that are only one strand and are no 5p/3p are determined as mature.
    counts_no_5p3p = counts_no_5p3p.drop('5p/3p', axis=1)
    rename_dict = dict(zip(libraries, libraries_mature))
    counts_no_5p3p = counts_no_5p3p.rename(columns=rename_dict)
    counts_no_5p3p['sum_FC_m'] = counts_no_5p3p[libraries_mature].sum(axis=1)
    print(counts_no_5p3p)
    mature_df = mature_df.append(counts_no_5p3p)

    counts_no_5p3p_filler = counts_no_5p3p.copy()
    rename_dict = dict(zip(libraries_mature, libraries_star))
    counts_no_5p3p_filler = counts_no_5p3p_filler.rename(columns=rename_dict)
    counts_no_5p3p_filler = counts_no_5p3p_filler.drop('sum_FC_m', axis=1)
    counts_no_5p3p_filler[libraries_star] = 0
    counts_no_5p3p_filler['sum_FC_s'] = 0
    counts_no_5p3p_filler['sum_FC_s > 100?'] = 0
    print(counts_no_5p3p_filler)
    star_df = star_df.append(counts_no_5p3p_filler)

    # Create index column for mirbase
    mirbase_intersections_table['index'] = mirbase_intersections_table['Description_mirbase'].str.split(';')
    mirbase_intersections_table['index'] = mirbase_intersections_table['index'].apply(lambda x : x[0])
    mirbase_intersections_table['index'] = mirbase_intersections_table['index'].str.replace('ID=MI', '')

    # Merge mirbase results and mirbase featurecounts results
    print(mirbase_intersections_table.loc[mirbase_intersections_table['index'] == '0007979'])
    print(mature_df.loc[mature_df['index'] == '0007979'])
    print(star_df.loc[star_df['index'] == '0007979'])
    mirbase_m_intersections_table = pd.merge(mirbase_intersections_table, mature_df, on='index', how='left')
    print(mirbase_m_intersections_table.loc[mirbase_m_intersections_table['index'] == '0007979'])
    mirbase_fc_intersections_table = pd.merge(mirbase_m_intersections_table, star_df, on='index', how='left')
    print(mirbase_fc_intersections_table.loc[mirbase_fc_intersections_table['index'] == '0007979'])
    mirbase_fc_intersections_table = mirbase_fc_intersections_table.drop('index', axis=1)
    print(mirbase_fc_intersections_table['Description_mirbase'].str.contains('0007979'))

    # filter by sum_fc_m < threshold
    mirbase_fc_intersections_table = mirbase_fc_intersections_table[mirbase_fc_intersections_table['sum_FC_m'] > sum_fc_thres]

    # Extract readcounts columns
    mirbase_fc_intersections_table['RC_m mirdeep'] = mirbase_fc_intersections_table["Description_mirdeep"].str.split(';', expand=True)[1]
    mirbase_fc_intersections_table['RC_s mirdeep'] = mirbase_fc_intersections_table["Description_mirdeep"].str.split(';', expand=True)[2]
    mirbase_fc_intersections_table['RC_m sRNAbench'] = mirbase_fc_intersections_table["Description_sRNAbench"].str.split(';', expand=True)[1]
    mirbase_fc_intersections_table['RC_s sRNAbench'] = mirbase_fc_intersections_table["Description_sRNAbench"].str.split(';', expand=True)[2]
    mirbase_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']] = mirbase_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']].fillna('0')
    mirbase_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']] = mirbase_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']].fillna('0')
    mirbase_fc_intersections_table['RC_m mirdeep'] = mirbase_fc_intersections_table['RC_m mirdeep'].str.replace('RC_m=', '').astype('int64')
    mirbase_fc_intersections_table['RC_s mirdeep'] = mirbase_fc_intersections_table['RC_s mirdeep'].str.replace('RC_s=', '').astype('int64')
    mirbase_fc_intersections_table['RC_m sRNAbench'] = mirbase_fc_intersections_table['RC_m sRNAbench'].str.replace('RC_m=', '').astype('int64')
    mirbase_fc_intersections_table['RC_s sRNAbench'] = mirbase_fc_intersections_table['RC_s sRNAbench'].str.replace('RC_s=', '').astype('int64')

    # Create diff columns
    mirbase_fc_intersections_table['Diff Sum_FC_m / RC_m mirdeep'] = mirbase_fc_intersections_table['sum_FC_m'] / mirbase_fc_intersections_table['RC_m mirdeep']
    mirbase_fc_intersections_table['Diff Sum_FC_m / RC_m sRNAbench'] = mirbase_fc_intersections_table['sum_FC_m'] / mirbase_fc_intersections_table['RC_m sRNAbench']
    mirbase_fc_intersections_table['Diff Sum_FC_s / RC_s mirdeep'] = mirbase_fc_intersections_table['sum_FC_s'] / mirbase_fc_intersections_table['RC_s mirdeep']
    mirbase_fc_intersections_table['Diff Sum_FC_s / RC_s sRNAbench'] = mirbase_fc_intersections_table['sum_FC_s'] / mirbase_fc_intersections_table['RC_s sRNAbench']

    # Normalize featurecounts to reads per million
    cast_dict = {k: 'int64' for k in libraries_mature + libraries_star}
    mirbase_fc_intersections_table = mirbase_fc_intersections_table.astype(cast_dict)
    total = mirbase_fc_intersections_table[libraries_mature + libraries_star].sum()
    mirbase_fc_intersections_table[mature_rpm + star_rpm] = round((mirbase_fc_intersections_table[libraries_mature + libraries_star] / total) * 1000000, 0)

    mirbase_fc_intersections_table['sum_FC_m_rpm'] = np.zeros(len(mirbase_fc_intersections_table))
    for library in mature_rpm:
        mirbase_fc_intersections_table['sum_FC_m_rpm'] += mirbase_fc_intersections_table[library]
    mirbase_fc_intersections_table['sum_FC_s_rpm'] = np.zeros(len(mirbase_fc_intersections_table))
    for library in star_rpm:
        mirbase_fc_intersections_table['sum_FC_s_rpm'] += mirbase_fc_intersections_table[library]
    mirbase_fc_intersections_table['mean_m_rpm'] = round(mirbase_fc_intersections_table[mature_rpm].mean(axis=1), 2)
    mirbase_fc_intersections_table['mean_s_rpm'] = round(mirbase_fc_intersections_table[star_rpm].mean(axis=1), 2)

# -----ADD SEQUENCES-----
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
mirdeep_blast_fc_intersections_table['mature_size'] = mirdeep_blast_fc_intersections_table['consensus mature sequence'].str.len()
mirdeep_blast_fc_intersections_table['star_size'] = mirdeep_blast_fc_intersections_table['consensus star sequence'].str.len()

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
sRNAbench_blast_fc_intersections_table['mature_size'] = np.where(sRNAbench_blast_fc_intersections_table['mature'] == '5p', sRNAbench_blast_fc_intersections_table['5pseq'].str.len(), sRNAbench_blast_fc_intersections_table['3pseq'].str.len())
sRNAbench_blast_fc_intersections_table['star_size'] = np.where(sRNAbench_blast_fc_intersections_table['mature'] == '5p', sRNAbench_blast_fc_intersections_table['3pseq'].str.len(), sRNAbench_blast_fc_intersections_table['5pseq'].str.len())


# -----REORDER COLUMNS:-----

if (species == 'Elegans') or (species == 'elegans'):
    elegans_columns_mirdeep = ['T/F_sRNAbench', 'Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase', 'T/F_mirbase',
                       'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb']
    elegans_columns_sRNAbench = ['T/F_mirdeep', 'Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase', 'T/F_mirbase',
                       'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb']
    mirbase_fc_intersections_table = mirbase_fc_intersections_table[['Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase',
                                                                     'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb',
                                                                     'Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep', 'T/F_mirdeep',
                                                                     'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench', 'T/F_sRNAbench'] +
                                                                    libraries_mature + ['sum_FC_m', 'RC_m mirdeep', 'RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench'] +
                                                                    libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s mirdeep', 'RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench'] +
                                                                    mature_rpm + ['sum_FC_m_rpm', 'mean_m_rpm'] + star_rpm + ['sum_FC_s_rpm', 'mean_s_rpm']
                                                                    ]
else:
    elegans_columns_mirdeep = []
    elegans_columns_sRNAbench = []


mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table[['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep',
                                                                             'query_accession','subject_accession', 'alignment_length', 'query_start', 'query_end', 'e_value',
                                                                             'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench'] + elegans_columns_mirdeep +
                                                                             libraries_mature + ['sum_FC_m', 'RC_m mirdeep', 'RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench'] +
                                                                             libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s mirdeep', 'RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench'] +
                                                                             mature_rpm + ['sum_FC_m_rpm', 'mean_m_rpm'] + star_rpm + ['sum_FC_s_rpm', 'mean_s_rpm'] +
                                                                             ['consensus mature sequence', 'consensus star sequence', 'consensus precursor sequence', 'mature', 'mature_size', 'star_size', 'loop_size']
                                                                            ]

sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table[['Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench',
                                                                             'query_accession','subject_accession', 'alignment_length', 'query_start', 'query_end', 'e_value',
                                                                             'Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep'] + elegans_columns_sRNAbench +
                                                                             libraries_mature + ['sum_FC_m', 'RC_m sRNAbench', 'RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep'] +
                                                                             libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s sRNAbench', 'RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep'] +
                                                                             mature_rpm + ['sum_FC_m_rpm', 'mean_m_rpm'] + star_rpm + ['sum_FC_s_rpm', 'mean_s_rpm'] +
                                                                             ['5pseq', '3pseq', 'hairpinSeq', 'mature', 'mature_size', 'star_size', 'loop_size']
                                                                                ]

# ------ADD TYPES------

if (species == 'Elegans') or (species == 'elegans'):
    # ------miRdeep:
    mirdeep_blast_fc_intersections_table['Type'] = np.zeros(len(mirdeep_blast_fc_intersections_table))
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 1
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 2
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 4
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 5
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 7
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 8


    # ------sRNAbench:
    sRNAbench_blast_fc_intersections_table['Type'] = np.zeros(len(sRNAbench_blast_fc_intersections_table))
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 1
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 2
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 4
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 5
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 7
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 8

    # ------mirbase:
    mirbase_fc_intersections_table['Type'] = np.zeros(len(mirbase_fc_intersections_table))
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 1
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 4
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 2
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 2
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 5
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 5
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 3
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 6

else:
    # ---miRdeep:
    mirdeep_blast_fc_intersections_table['Type'] = np.where(mirdeep_blast_fc_intersections_table['Description_sRNAbench'] != '.', 1, 2)

    # ---sRNAbench:
    sRNAbench_blast_fc_intersections_table['Type'] = np.where(sRNAbench_blast_fc_intersections_table['Description_mirdeep'] != '.', 1, 3)


# ------CREATE ALL CANDIDATES SHEET:------
unified = mirdeep_blast_fc_intersections_table.copy()

# Defining 5p/3p sequences in mirdeep
unified['5pseq'] = np.where(unified['mature'] == '5p', unified['consensus mature sequence'], unified['consensus star sequence'])
unified['3pseq'] = np.where(unified['mature'] == '3p', unified['consensus mature sequence'], unified['consensus star sequence'])

unified = unified.drop(['Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench', 'consensus mature sequence', 'consensus star sequence'], axis=1)
if (species == 'Elegans') or (species == 'elegans'):
    only_sRNAbench = sRNAbench_blast_fc_intersections_table[sRNAbench_blast_fc_intersections_table['Type'].isin([2, 5, 8])]
    only_mirbase = mirbase_fc_intersections_table[mirbase_fc_intersections_table['Type'].isin([3, 6])]
   # unified = unified.drop(['T/F_sRNAbench', 'T/F_mirbase', 'T/F_mirgenedb'], axis=1)
    #only_sRNAbench = only_sRNAbench.drop(['T/F_mirdeep'], axis=1)
else:
    only_sRNAbench = sRNAbench_blast_fc_intersections_table[sRNAbench_blast_fc_intersections_table['Type'] == 3]


only_sRNAbench = only_sRNAbench.drop(['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep'], axis=1)
unified = unified.rename(columns={'Chr_mirdeep':'Chr', 'Start_mirdeep':'Start', 'End_mirdeep':'End', 'Strand_mirdeep':'Strand', 'Description_mirdeep':'Description', 'consensus precursor sequence':'hairpinSeq'})
only_sRNAbench = only_sRNAbench.rename(columns={'Chr_sRNAbench':'Chr', 'Start_sRNAbench':'Start', 'End_sRNAbench':'End', 'Strand_sRNAbench':'Strand', 'Description_sRNAbench':'Description'})
if (species == 'Elegans') or (species == 'elegans'):
    only_mirbase = only_mirbase.drop(['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep', 'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench',], axis=1)
    only_mirbase = only_mirbase.rename(columns={'Chr_mirbase':'Chr', 'Start_mirbase':'Start', 'End_mirbase':'End', 'Strand_mirbase':'Strand', 'Description_mirbase':'Description'})
    unified = pd.concat([unified, only_sRNAbench, only_mirbase], axis=0, ignore_index=True)
    unified = unified.drop(['T/F_mirdeep', 'T/F_sRNAbench', 'T/F_mirbase', 'T/F_mirgenedb'], axis=1)
else:
    unified = pd.concat([unified, only_sRNAbench], axis=0, ignore_index=True)
# reorder sequences columns
columns = list(unified.columns)
i1, i2, i3 = columns.index('mature'), columns.index('5pseq'), columns.index('3pseq')
columns[i1], columns[i2], columns[i3] = columns[i2], columns[i3], columns[i1]
unified = unified[columns]


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


# -----STATISTICS-----
# ---Creating a families by type pivot table


def families_by_type(df):
    families_by_type = pd.pivot_table(df, values='Description', index='Family', columns='Type', aggfunc='count')
    families_by_type.fillna(0, inplace=True)
    families_by_type = families_by_type.drop("UNKNOWN", axis=0)
    families_by_type.sort_values('Family', inplace=True)
    families_by_type.plot(kind='bar', figsize=(14, 10), stacked=True, title="{} counts of known families in each type".format(species))
    plt.ylabel("Counts")
    # plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
    # plt.show()
    plt.savefig("{}_family_counts_per_type.png".format(species))


families_by_type(unified)

# ---Creating unknown seeds by type pivot table


def unknown_families_by_type(df):
    filtered_df = df[df['Family'] == "UNKNOWN"]
    pivot = pd.pivot_table(filtered_df, values='Description', index='Seed', columns='Type', aggfunc='count')
    pivot.fillna(0, inplace=True)
    pivot['sum'] = pivot.sum(axis=1)
    pivot_groups = pivot[pivot['sum'] > 1]
    pivot_solos = pivot[pivot['sum'] == 1]

    pivot_groups = pivot_groups.drop('sum', axis=1)
    pivot_groups.sort_values('Seed', inplace=True)
    pivot_groups.plot(kind='bar', figsize=(14, 10), stacked=True,
                          title="{} counts of unknown families in each type".format(species))
    plt.ylabel("Counts")
    # plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
    # plt.show()
    plt.savefig("{}_unknown_family_counts_per_type.png".format(species))
    plt.clf()

    pivot_solos = pivot_solos.drop('sum', axis=1)
    pivot_solos.sum().plot(kind='pie', y='Type', autopct='%1.0f%%',
                           title="{} unique seed candidates by type".format(species))
    plt.ylabel("Counts")
    # plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
    # plt.show()
    plt.savefig("{}_unique_seed_candidates_by_type.png".format(species))



unknown_families_by_type(unified)

# -----SAVE TO EXCEL-----

writer = pd.ExcelWriter('intersections_table_{}.xlsx'.format(species))
mirdeep_blast_fc_intersections_table.to_excel(writer, sheet_name='miRdeep')
sRNAbench_blast_fc_intersections_table.to_excel(writer, sheet_name='sRNAbench')
if (species == 'Elegans') or (species == 'elegans'):
    mirbase_fc_intersections_table.to_excel(writer, sheet_name='mirbase')
unified.to_excel(writer, sheet_name='all_candidates')
blast_mirdeep_orig.to_excel(writer, sheet_name='blast_miRdeep')
blast_sRNAbench_orig.to_excel(writer, sheet_name='blast_sRNAbench')
writer.save()

