import numpy as np
import pandas as pd

"""
Create intersections table
Add BLAST results
Add feature counts
Add sequences
Add types
"""

# -----mirdeep intersections table:-----

mirdeep_sRNAbench = pd.read_csv('miRdeep_sRNAbench_intersect.bed', sep='\t', names=['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2', 'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_sRNAbench', '.4', 'pre_miRNA2', 'Start_sRNAbench', 'End_sRNAbench', '.5', 'Strand_sRNAbench', '.6', 'Description_sRNAbench'])
mirdeep_mirbase = pd.read_csv('miRdeep_miRBase_intersect.bed', sep='\t', names=['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2', 'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_mirbase', '.4', 'pre_miRNA2', 'Start_mirbase', 'End_mirbase', '.5', 'Strand_mirbase', '.6', 'Description_mirbase'])
mirdeep_mirgenedb = pd.read_csv('miRdeep_miRGeneDB_intersect.bed', sep='\t', names=['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2', 'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_mirgenedb', '.4', 'pre_miRNA2', 'Start_mirgenedb', 'End_mirgenedb', '.5', 'Strand_mirgenedb', '.6', 'Description_mirgenedb'])

mirdeep_sRNAbench = mirdeep_sRNAbench.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
mirdeep_mirbase = mirdeep_mirbase.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
mirdeep_mirgenedb = mirdeep_mirgenedb.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)

mirdeep_sRNAbench['T/F_sRNAbench'] = (mirdeep_sRNAbench['Description_sRNAbench'] != '.').astype(int) # Used for classifying types

mirdeep_sRNAbench_mirbase = pd.merge(mirdeep_sRNAbench, mirdeep_mirbase.iloc[:, 4:10], on='Description_mirdeep', how='left')
mirdeep_sRNAbench_mirbase['T/F_mirbase'] = (mirdeep_sRNAbench_mirbase['Description_mirbase'] != '.').astype(int) # Used for classifying types
mirdeep_intersections_table = pd.merge(mirdeep_sRNAbench_mirbase, mirdeep_mirgenedb.iloc[:, 4:10], on='Description_mirdeep', how='left')
mirdeep_intersections_table['T/F_mirgenedb'] = (mirdeep_intersections_table['Description_mirgenedb'] != '.').astype(int) # Used for classifying types
mirdeep_intersections_table.to_csv("mirdeep_intersections_table", sep='\t')

# -----sRNAbench intersections table:-----

sRNAbench_mirdeep = pd.read_csv('sRNAbench_miRdeep_intersect.bed', sep='\t', names=['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench', '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirdeep', '.4', 'pre_miRNA2', 'Start_mirdeep', 'End_mirdeep', '.5', 'Strand_mirdeep', '.6', 'Description_mirdeep'])
sRNAbench_mirbase = pd.read_csv('sRNAbench_miRBase_intersect.bed', sep='\t', names=['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench', '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirbase', '.4', 'pre_miRNA2', 'Start_mirbase', 'End_mirbase', '.5', 'Strand_mirbase', '.6', 'Description_mirbase'])
sRNAbench_mirgenedb = pd.read_csv('sRNAbench_miRGeneDB_intersect.bed', sep='\t', names=['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench', '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirgenedb', '.4', 'pre_miRNA2', 'Start_mirgenedb', 'End_mirgenedb', '.5', 'Strand_mirgenedb', '.6', 'Description_mirgenedb'])

sRNAbench_mirdeep = sRNAbench_mirdeep.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
sRNAbench_mirbase = sRNAbench_mirbase.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
sRNAbench_mirgenedb = sRNAbench_mirgenedb.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)

sRNAbench_mirdeep['T/F_mirdeep'] = (sRNAbench_mirdeep['Description_mirdeep'] != '.').astype(int) # Used for classifying types

sRNAbench_mirdeep_mirbase = pd.merge(sRNAbench_mirdeep, sRNAbench_mirbase.iloc[:, 4:10], on='Description_sRNAbench', how='left')
sRNAbench_mirdeep_mirbase['T/F_mirbase'] = (sRNAbench_mirdeep_mirbase['Description_mirbase'] != '.').astype(int) # Used for classifying types
sRNAbench_intersections_table = pd.merge(sRNAbench_mirdeep_mirbase, sRNAbench_mirgenedb.iloc[:, 4:10], on='Description_sRNAbench', how='left')
sRNAbench_intersections_table['T/F_mirgenedb'] = (sRNAbench_intersections_table['Description_mirgenedb'] != '.').astype(int) # Used for classifying types
sRNAbench_intersections_table.to_csv("sRNAbench_intersections_table", sep='\t')

# -----mirbase intersections table:-----

mirbase_mirgenedb = pd.read_csv('miRBase_miRGeneDB_intersect.bed', sep='\t', names=['Chr_mirbase', '.1', 'pre_miRNA1', 'Start_mirbase', 'End_mirbase', '.2', 'Strand_mirbase', '.3', 'Description_mirbase', 'Chr_mirgenedb', '.4', 'pre_miRNA2', 'Start_mirgenedb', 'End_mirgenedb', '.5', 'Strand_mirgenedb', '.6', 'Description_mirgenedb'])
mirbase_mirdeep = pd.read_csv('miRBase_miRdeep_intersect.bed', sep='\t', names=['Chr_mirbase', '.1', 'pre_miRNA1', 'Start_mirbase', 'End_mirbase', '.2', 'Strand_mirbase', '.3', 'Description_mirbase', 'Chr_mirdeep', '.4', 'pre_miRNA2', 'Start_mirdeep', 'End_mirdeep', '.5', 'Strand_mirdeep', '.6', 'Description_mirdeep'])
mirbase_sRNAbench = pd.read_csv('miRBase_sRNAbench_intersect.bed', sep='\t', names=['Chr_mirbase', '.1', 'pre_miRNA1', 'Start_mirbase', 'End_mirbase', '.2', 'Strand_mirbase', '.3', 'Description_mirbase', 'Chr_sRNAbench', '.4', 'pre_miRNA2', 'Start_sRNAbench', 'End_sRNAbench', '.5', 'Strand_sRNAbench', '.6', 'Description_sRNAbench'])

mirbase_mirgenedb = mirbase_mirgenedb.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
mirbase_mirdeep = mirbase_mirdeep.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
mirbase_sRNAbench = mirbase_sRNAbench.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)

mirbase_mirgenedb['T/F_mirgenedb'] = (mirbase_mirgenedb['Description_mirgenedb'] != '.').astype(int) # Used for classifying types

mirbase_mirgenedb_mirdeep = pd.merge(mirbase_mirgenedb, mirbase_mirdeep.iloc[:, 4:10], on='Description_mirbase', how='left')
mirbase_mirgenedb_mirdeep['T/F_mirdeep'] = (mirbase_mirgenedb_mirdeep['Description_mirdeep'] != '.').astype(int) # Used for classifying types
mirbase_intersections_table = pd.merge(mirbase_mirgenedb_mirdeep, mirbase_sRNAbench.iloc[:, 4:10], on='Description_mirbase', how='left')
mirbase_intersections_table['T/F_sRNAbench'] = (mirbase_intersections_table['Description_sRNAbench'] != '.').astype(int) # Used for classifying types
mirbase_intersections_table.to_csv("mirbase_intersections_table", sep='\t')

# -----Add blast results:-----
# -----miRdeep:
blast_mirdeep_orig = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/queries/Elegans/miRdeep_blastn_compact', sep='\t', names=['query_accession', 'subject_accession', '%_identical_matches', 'alignment_length', 'mismatches', 'gap_openings', 'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bitscore'])
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

# -----sRNAbench:

blast_sRNAbench_orig = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/queries/Elegans/sRNAbench_blastn_compact', sep='\t', names=['query_accession', 'subject_accession', '%_identical_matches', 'alignment_length', 'mismatches', 'gap_openings', 'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bitscore'])
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
# -----miRDeep:
featurecounts_mirdeep = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_miRdeep_counts.txt', sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'CE57', 'CE58', 'CE59', 'CE60', 'CE61', 'CE62', 'CE63', 'CE69', 'CE78', 'CE79', 'CE80', 'CE81'])
featurecounts_mirdeep = featurecounts_mirdeep.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
featurecounts_mirdeep = featurecounts_mirdeep.iloc[2:] # Drop the first 2 rows, which is readme info from featurecounts and not data

# Create mature/star and index column
featurecounts_mirdeep['index'] = featurecounts_mirdeep['Geneid'].str.split('|')
featurecounts_mirdeep['index'] = featurecounts_mirdeep['index'].apply(lambda x : x[4])

featurecounts_mirdeep['mature/star'] = featurecounts_mirdeep['Geneid'].str.split('|')
featurecounts_mirdeep['mature/star'] = featurecounts_mirdeep['mature/star'].apply(lambda x : x[2])

featurecounts_mirdeep = featurecounts_mirdeep.drop('Geneid', axis=1)

# Separate df into mature and star
mature_counts = featurecounts_mirdeep[featurecounts_mirdeep['mature/star'] == 'm']
mature_counts = mature_counts.rename(columns={'CE57':'CE57_m', 'CE58':'CE58_m', 'CE59':'CE59_m', 'CE60':'CE60_m', 'CE61':'CE61_m', 'CE62':'CE62_m', 'CE63':'CE63_m', 'CE69':'CE69_m', 'CE78':'CE78_m', 'CE79':'CE79_m', 'CE80':'CE80_m', 'CE81':'CE81_m'})
mature_counts = mature_counts.drop('mature/star', axis=1)

star_counts = featurecounts_mirdeep[featurecounts_mirdeep['mature/star'] == 's']
star_counts = star_counts.rename(columns={'CE57':'CE57_s', 'CE58':'CE58_s', 'CE59':'CE59_s', 'CE60':'CE60_s', 'CE61':'CE61_s', 'CE62':'CE62_s', 'CE63':'CE63_s', 'CE69':'CE69_s', 'CE78':'CE78_s', 'CE79':'CE79_s', 'CE80':'CE80_s', 'CE81':'CE81_s'})
star_counts = star_counts.drop('mature/star', axis=1)

# Merge mirdeep & blast results and featurecounts results
mirdeep_blast_m_intersections_table = pd.merge(mirdeep_blast_intersections_table, mature_counts, on='index', how='left')
mirdeep_blast_fc_intersections_table = pd.merge(mirdeep_blast_m_intersections_table, star_counts, on='index', how='left')
# mirdeep_blast_fc_intersections_table = pd.merge(mirdeep_blast_fc_intersections_table, mature_counts_mb, on='index', how='left')
# mirdeep_blast_fc_intersections_table = pd.merge(mirdeep_blast_fc_intersections_table, star_counts_mb, on='index', how='left')
mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table.drop('index', axis=1)

# -----sRNAbench:
featurecounts_sRNAbench = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_sRNAbench_counts.txt', sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'CE57', 'CE58', 'CE59', 'CE60', 'CE61', 'CE62', 'CE63', 'CE69', 'CE78', 'CE79', 'CE80', 'CE81'])
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

# Separate df into mature and star
mature_counts = featurecounts_sRNAbench[featurecounts_sRNAbench['mature/star'] == 'm']
mature_counts = mature_counts.rename(columns={'CE57':'CE57_m', 'CE58':'CE58_m', 'CE59':'CE59_m', 'CE60':'CE60_m', 'CE61':'CE61_m', 'CE62':'CE62_m', 'CE63':'CE63_m', 'CE69':'CE69_m', 'CE78':'CE78_m', 'CE79':'CE79_m', 'CE80':'CE80_m', 'CE81':'CE81_m'})
mature_counts = mature_counts.drop('mature/star', axis=1)

star_counts = featurecounts_sRNAbench[featurecounts_sRNAbench['mature/star'] == 's']
star_counts = star_counts.rename(columns={'CE57':'CE57_s', 'CE58':'CE58_s', 'CE59':'CE59_s', 'CE60':'CE60_s', 'CE61':'CE61_s', 'CE62':'CE62_s', 'CE63':'CE63_s', 'CE69':'CE69_s', 'CE78':'CE78_s', 'CE79':'CE79_s', 'CE80':'CE80_s', 'CE81':'CE81_s'})
star_counts = star_counts.drop(['mature/star', 'mature'], axis=1)

# Merge sRNAbench & blast results and featurecounts results
sRNAbench_blast_m_intersections_table = pd.merge(sRNAbench_blast_intersections_table, mature_counts, on='index', how='left')
sRNAbench_blast_fc_intersections_table = pd.merge(sRNAbench_blast_m_intersections_table, star_counts, on='index', how='left')
# sRNAbench_blast_fc_intersections_table = pd.merge(sRNAbench_blast_fc_intersections_table, mature_counts_mb, on='index', how='left')
# sRNAbench_blast_fc_intersections_table = pd.merge(sRNAbench_blast_fc_intersections_table, star_counts_mb, on='index', how='left')
sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table.drop('index', axis=1)

# -----miRBase:
featurecounts_mirbase = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Elegans/counts_sep/miRNA_mirbase_counts.txt', sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'mb_CE57', 'mb_CE58', 'mb_CE59', 'mb_CE60', 'mb_CE61', 'mb_CE62', 'mb_CE63', 'mb_CE69', 'mb_CE78', 'mb_CE79', 'mb_CE80', 'mb_CE81'])
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

# Separate df into 5p and 3p
counts_mb_5p = featurecounts_mirbase[featurecounts_mirbase['5p/3p'] == '5p']
counts_mb_5p = counts_mb_5p.rename(columns={'mb_CE57':'mb_CE57_5p', 'mb_CE58':'mb_CE58_5p', 'mb_CE59':'mb_CE59_5p', 'mb_CE60':'mb_CE60_5p', 'mb_CE61':'mb_CE61_5p', 'mb_CE62':'mb_CE62_5p', 'mb_CE63':'mb_CE63_5p', 'mb_CE69':'mb_CE69_5p', 'mb_CE78':'mb_CE78_5p', 'mb_CE79':'mb_CE79_5p', 'mb_CE80':'mb_CE80_5p', 'mb_CE81':'mb_CE81_5p'})
counts_mb_5p = counts_mb_5p.drop('5p/3p', axis=1)

counts_mb_3p = featurecounts_mirbase[featurecounts_mirbase['5p/3p'] == '3p']
counts_mb_3p = counts_mb_3p.rename(columns={'mb_CE57':'mb_CE57_3p', 'mb_CE58':'mb_CE58_3p', 'mb_CE59':'mb_CE59_3p', 'mb_CE60':'mb_CE60_3p', 'mb_CE61':'mb_CE61_3p', 'mb_CE62':'mb_CE62_3p', 'mb_CE63':'mb_CE63_3p', 'mb_CE69':'mb_CE69_3p', 'mb_CE78':'mb_CE78_3p', 'mb_CE79':'mb_CE79_3p', 'mb_CE80':'mb_CE80_3p', 'mb_CE81':'mb_CE81_3p'})
counts_mb_3p = counts_mb_3p.drop('5p/3p', axis=1)

counts_no_5p3p = featurecounts_mirbase[(featurecounts_mirbase['5p/3p'] != '3p') & (featurecounts_mirbase['5p/3p'] != '5p')]

# sum 5p and 3p to determine mature/star
numeric_counts_mb_5p = counts_mb_5p[['mb_CE57_5p', 'mb_CE58_5p', 'mb_CE59_5p', 'mb_CE60_5p', 'mb_CE61_5p', 'mb_CE62_5p', 'mb_CE63_5p', 'mb_CE69_5p', 'mb_CE78_5p', 'mb_CE79_5p', 'mb_CE80_5p', 'mb_CE81_5p']].astype(int)
counts_mb_5p['sum'] = numeric_counts_mb_5p.sum(axis=1)
numeric_counts_mb_3p = counts_mb_3p[['mb_CE57_3p', 'mb_CE58_3p', 'mb_CE59_3p', 'mb_CE60_3p', 'mb_CE61_3p', 'mb_CE62_3p', 'mb_CE63_3p', 'mb_CE69_3p', 'mb_CE78_3p', 'mb_CE79_3p', 'mb_CE80_3p', 'mb_CE81_3p']].astype(int)
counts_mb_3p['sum'] = numeric_counts_mb_3p.sum(axis=1)

# Create empty mature df and star df, iterate the rows of 5p and 3p and add to the relavant.
mature_df = pd.DataFrame(columns=['mb_CE57_m', 'mb_CE58_m', 'mb_CE59_m', 'mb_CE60_m', 'mb_CE61_m', 'mb_CE62_m', 'mb_CE63_m', 'mb_CE69_m', 'mb_CE78_m', 'mb_CE79_m', 'mb_CE80_m', 'mb_CE81_m'])
star_df = pd.DataFrame(columns=['mb_CE57_s', 'mb_CE58_s', 'mb_CE59_s', 'mb_CE60_s', 'mb_CE61_s', 'mb_CE62_s', 'mb_CE63_s', 'mb_CE69_s', 'mb_CE78_s', 'mb_CE79_s', 'mb_CE80_s', 'mb_CE81_s'])
for i in range(0, len(counts_mb_5p)):
    row_5p = counts_mb_5p.iloc[i]
    row_3p = counts_mb_3p.iloc[i]
    if row_5p['sum'] > row_3p['sum']: # if 5p is the mature
        row_5p = row_5p.rename({'mb_CE57_5p': 'mb_CE57_m', 'mb_CE58_5p': 'mb_CE58_m', 'mb_CE59_5p': 'mb_CE59_m',
                                'mb_CE60_5p': 'mb_CE60_m', 'mb_CE61_5p': 'mb_CE61_m', 'mb_CE62_5p': 'mb_CE62_m',
                                'mb_CE63_5p': 'mb_CE63_m', 'mb_CE69_5p': 'mb_CE69_m', 'mb_CE78_5p': 'mb_CE78_m',
                                'mb_CE79_5p': 'mb_CE79_m', 'mb_CE80_5p': 'mb_CE80_m',
                                'mb_CE81_5p': 'mb_CE81_m'})
        row_3p = row_3p.rename({'mb_CE57_3p': 'mb_CE57_s', 'mb_CE58_3p': 'mb_CE58_s', 'mb_CE59_3p': 'mb_CE59_s',
                                'mb_CE60_3p': 'mb_CE60_s', 'mb_CE61_3p': 'mb_CE61_s', 'mb_CE62_3p': 'mb_CE62_s',
                                'mb_CE63_3p': 'mb_CE63_s', 'mb_CE69_3p': 'mb_CE69_s', 'mb_CE78_3p': 'mb_CE78_s',
                                'mb_CE79_3p': 'mb_CE79_s', 'mb_CE80_3p': 'mb_CE80_s',
                                'mb_CE81_3p': 'mb_CE81_s'})
        mature_df = mature_df.append(row_5p)
        star_df = star_df.append(row_3p)
    else: # else 3p is the mature
        row_5p = row_5p.rename({'mb_CE57_5p': 'mb_CE57_s', 'mb_CE58_5p': 'mb_CE58_s', 'mb_CE59_5p': 'mb_CE59_s',
                                'mb_CE60_5p': 'mb_CE60_s', 'mb_CE61_5p': 'mb_CE61_s', 'mb_CE62_5p': 'mb_CE62_s',
                                'mb_CE63_5p': 'mb_CE63_s', 'mb_CE69_5p': 'mb_CE69_s', 'mb_CE78_5p': 'mb_CE78_s',
                                'mb_CE79_5p': 'mb_CE79_s', 'mb_CE80_5p': 'mb_CE80_s',
                                'mb_CE81_5p': 'mb_CE81_s'})
        row_3p = row_3p.rename({'mb_CE57_3p': 'mb_CE57_m', 'mb_CE58_3p': 'mb_CE58_m', 'mb_CE59_3p': 'mb_CE59_m',
                                'mb_CE60_3p': 'mb_CE60_m', 'mb_CE61_3p': 'mb_CE61_m', 'mb_CE62_3p': 'mb_CE62_m',
                                'mb_CE63_3p': 'mb_CE63_m', 'mb_CE69_3p': 'mb_CE69_m', 'mb_CE78_3p': 'mb_CE78_m',
                                'mb_CE79_3p': 'mb_CE79_m', 'mb_CE80_3p': 'mb_CE80_m',
                                'mb_CE81_3p': 'mb_CE81_m'})
        mature_df = mature_df.append(row_3p)
        star_df = star_df.append(row_5p)

mature_df = mature_df.drop('sum', axis=1)
star_df = star_df.drop('sum', axis=1)

# Those that are only one strand and are no 5p/3p are determined as mature.
counts_no_5p3p = counts_no_5p3p.drop('5p/3p', axis=1)
counts_no_5p3p = counts_no_5p3p.rename(columns={'mb_CE57': 'mb_CE57_m', 'mb_CE58': 'mb_CE58_m', 'mb_CE59': 'mb_CE59_m',
                                                'mb_CE60': 'mb_CE60_m', 'mb_CE61': 'mb_CE61_m', 'mb_CE62': 'mb_CE62_m',
                                                'mb_CE63': 'mb_CE63_m', 'mb_CE69': 'mb_CE69_m', 'mb_CE78': 'mb_CE78_m',
                                                'mb_CE79': 'mb_CE79_m', 'mb_CE80': 'mb_CE80_m', 'mb_CE81': 'mb_CE81_m'})
mature_df = mature_df.append(counts_no_5p3p)

# Create index column for mirbase
mirbase_intersections_table['index'] = mirbase_intersections_table['Description_mirbase'].str.split(';')
mirbase_intersections_table['index'] = mirbase_intersections_table['index'].apply(lambda x : x[0])
mirbase_intersections_table['index'] = mirbase_intersections_table['index'].str.replace('ID=MI', '')

# Merge mirbase results and mirbase featurecounts results
mirbase_m_intersections_table = pd.merge(mirbase_intersections_table, mature_df, on='index', how='left')
mirbase_fc_intersections_table = pd.merge(mirbase_m_intersections_table, star_df, on='index', how='left')
mirbase_fc_intersections_table = mirbase_fc_intersections_table.drop('index', axis=1)

# -----Reorder columns:-----
mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table[['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep',
                                                                             'query_accession','subject_accession', 'alignment_length', 'query_start', 'query_end', 'strand', 'e_value',
                                                                             'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench', 'T/F_sRNAbench',
                                                                             'Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase', 'T/F_mirbase',
                                                                             # 'mb_CE57_m', 'mb_CE58_m', 'mb_CE59_m', 'mb_CE60_m', 'mb_CE61_m', 'mb_CE62_m', 'mb_CE63_m', 'mb_CE69_m', 'mb_CE78_m', 'mb_CE79_m', 'mb_CE80_m', 'mb_CE81_m',
                                                                             # 'mb_CE57_s', 'mb_CE58_s', 'mb_CE59_s', 'mb_CE60_s', 'mb_CE61_s', 'mb_CE62_s', 'mb_CE63_s', 'mb_CE69_s', 'mb_CE78_s', 'mb_CE79_s', 'mb_CE80_s', 'mb_CE81_s',
                                                                             'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb',
                                                                             'CE57_m', 'CE58_m','CE59_m','CE60_m','CE61_m',	'CE62_m', 'CE63_m', 'CE69_m', 'CE78_m', 'CE79_m', 'CE80_m', 'CE81_m',
                                                                             'CE57_s', 'CE58_s', 'CE59_s', 'CE60_s', 'CE61_s', 'CE62_s', 'CE63_s', 'CE69_s', 'CE78_s', 'CE79_s', 'CE80_s', 'CE81_s'
                                                                             ]]

sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table[['Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench',
                                                                             'query_accession','subject_accession', 'alignment_length', 'query_start', 'query_end', 'strand', 'e_value',
                                                                             'Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep', 'T/F_mirdeep',
                                                                             'Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase', 'T/F_mirbase',
                                                                             # 'mb_CE57_m', 'mb_CE58_m', 'mb_CE59_m', 'mb_CE60_m', 'mb_CE61_m', 'mb_CE62_m', 'mb_CE63_m', 'mb_CE69_m', 'mb_CE78_m', 'mb_CE79_m', 'mb_CE80_m', 'mb_CE81_m',
                                                                             # 'mb_CE57_s', 'mb_CE58_s', 'mb_CE59_s', 'mb_CE60_s', 'mb_CE61_s', 'mb_CE62_s', 'mb_CE63_s', 'mb_CE69_s', 'mb_CE78_s', 'mb_CE79_s', 'mb_CE80_s', 'mb_CE81_s',
                                                                             'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb',
                                                                             'CE57_m', 'CE58_m','CE59_m','CE60_m','CE61_m',	'CE62_m', 'CE63_m', 'CE69_m', 'CE78_m', 'CE79_m', 'CE80_m', 'CE81_m',
                                                                             'CE57_s', 'CE58_s', 'CE59_s', 'CE60_s', 'CE61_s', 'CE62_s', 'CE63_s', 'CE69_s', 'CE78_s', 'CE79_s', 'CE80_s', 'CE81_s', 'mature'
                                                                             ]]

# -----Add Sequences:-----
# -----miRdeep:
remaining1_mirdeep = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Elegans/mirdeep_out/remaining_file_1.csv', sep='\t')
remaining2_mirdeep = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Elegans/mirdeep_out/remaining_file_2.csv', sep='\t')
remaining_mirdeep = pd.concat([remaining1_mirdeep, remaining2_mirdeep], ignore_index=True)
mirdeep_blast_fc_intersections_table['consensus mature sequence'] = remaining_mirdeep['consensus mature sequence']
mirdeep_blast_fc_intersections_table['consensus star sequence'] = remaining_mirdeep['consensus star sequence']
mirdeep_blast_fc_intersections_table['consensus precursor sequence'] = remaining_mirdeep['consensus precursor sequence']
#-----sRNAbench:
remaining_sRNAbench = pd.read_csv('/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Elegans/sRNAbench_remaining.csv', sep='\t')
sRNAbench_blast_fc_intersections_table['5pseq'] = remaining_sRNAbench['5pseq']
sRNAbench_blast_fc_intersections_table['3pseq'] = remaining_sRNAbench['3pseq']
sRNAbench_blast_fc_intersections_table['hairpinSeq'] = remaining_sRNAbench['hairpinSeq']

sRNAbench_blast_fc_intersections_table['Seed'] = np.where(sRNAbench_blast_fc_intersections_table["mature"] == '5p', sRNAbench_blast_fc_intersections_table["5pseq"].str[1:8], sRNAbench_blast_fc_intersections_table["3pseq"].str[1:8])

# ------Add Types------
# ------miRdeep:
mirdeep_blast_fc_intersections_table['Type'] = np.zeros(len(mirdeep_blast_fc_intersections_table))
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 1
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 2
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 3
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 4
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 5
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 6
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 7
mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 8

# ------sRNAbench:
sRNAbench_blast_fc_intersections_table['Type'] = np.zeros(len(sRNAbench_blast_fc_intersections_table))
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 1
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 2
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 3
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 4
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 5
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 6
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 7
sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 8

# ------mirbase:
mirbase_fc_intersections_table['Type'] = np.zeros(len(mirbase_fc_intersections_table))
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 1
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 2
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 3
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 4
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 5
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 6
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 7
mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 8

# -----Create unified sheet:-----
unified = mirbase_fc_intersections_table.copy()
unified = unified.drop([], axis=1)

# -----Save each intersections table as a sheet in one excel file:-----

writer = pd.ExcelWriter('intersections_table_elegans.xlsx')
mirdeep_blast_fc_intersections_table.to_excel(writer, sheet_name='miRdeep')
sRNAbench_blast_fc_intersections_table.to_excel(writer, sheet_name='sRNAbench')
mirbase_fc_intersections_table.to_excel(writer, sheet_name='mirbase')
blast_mirdeep_orig.to_excel(writer, sheet_name='blast_miRdeep')
blast_sRNAbench_orig.to_excel(writer, sheet_name='blast_sRNAbench')
writer.save()

