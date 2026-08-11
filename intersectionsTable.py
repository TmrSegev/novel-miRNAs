import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import gffpandas.gffpandas as gffpd

from pipeline_config import get_species_config, load_seed_table, species_uses_blast, species_uses_mirbase
from pipeline_blast import load_blast_results, merge_blast

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


def normalize_description(series):
    """Coerce bedtools/GFF description fields to string; nulls / empty → '.'."""
    return (
        series.fillna(".")
        .replace("", ".")
        .astype(str)
        .replace({"nan": ".", "None": ".", "<NA>": "."})
    )


def extract_description_field(series, field_index):
    """Split description on ';'; return field_index, or NaN if that field is absent.

    Needed when -loj leaves Description as '.' / empty for non-overlaps: pandas may
    type an all-null column as float64 (breaks .str), and expand may have < field_index+1 cols.
    """
    parts = normalize_description(series).str.split(";", expand=True)
    if parts.shape[1] <= field_index:
        return pd.Series(np.nan, index=series.index)
    return parts[field_index]


def load_cross_intersections(path, columns, a_tool, b_tool):
    """Load a full A+B bedtools intersection and reject incomplete output."""
    expected_fields = len(columns)
    observed_fields = set()
    with open(path) as handle:
        for line in handle:
            if line.strip():
                observed_fields.add(len(line.rstrip("\n").split("\t")))

    if not observed_fields:
        raise ValueError("Cross-tool intersection file is empty: {}".format(path))
    if observed_fields != {expected_fields}:
        raise ValueError(
            "{} -> {} intersection must have {} tab-separated fields per row, "
            "but found {} in {}. Re-run Phase 9 with "
            "'bedtools intersect -wa -wb -loj -s -f 0.6'.".format(
                a_tool, b_tool, expected_fields, sorted(observed_fields), path
            )
        )

    table = pd.read_csv(path, sep="\t", names=columns)
    table["Description_{}".format(a_tool)] = normalize_description(
        table["Description_{}".format(a_tool)]
    )
    table["Description_{}".format(b_tool)] = normalize_description(
        table["Description_{}".format(b_tool)]
    )
    return table


def attach_remaining_columns(table, remaining, description_col, column_map, tool):
    """Attach candidate metadata using the GFF index field, never row position."""
    index_field = extract_description_field(table[description_col], 3)
    candidate_index = pd.to_numeric(
        index_field.astype(str).str.replace("index=", ""), errors="coerce"
    )
    if candidate_index.isna().any():
        raise ValueError(
            "{} primary descriptions contain missing/invalid index fields.".format(tool)
        )

    candidate_index = candidate_index.astype("int64")
    expected = set(range(len(remaining)))
    observed = set(candidate_index.tolist())
    if observed != expected:
        missing = sorted(expected - observed)
        unexpected = sorted(observed - expected)
        raise ValueError(
            "{} intersection candidates do not match the remaining-candidates file: "
            "{} rows in remaining, {} unique intersection indices; "
            "missing indices {}{}; unexpected indices {}{}. "
            "Re-run Phase 5 and then Phase 9 with -wa -wb -loj.".format(
                tool,
                len(remaining),
                len(observed),
                missing[:10],
                "..." if len(missing) > 10 else "",
                unexpected[:10],
                "..." if len(unexpected) > 10 else "",
            )
        )

    remaining_by_index = remaining.reset_index(drop=True)
    for output_col, source_col in column_map.items():
        table[output_col] = candidate_index.map(remaining_by_index[source_col])
    return table

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
featurecounts_pre_mirdeep_path = None
featurecounts_pre_sRNAbench_path = None
featurecounts_mirbase_path = None
remaining_mirdeep_path = None
remaining1_mirdeep_path = None
remaining2_mirdeep_path = None
remaining_sRNAbench_path = None
mirbase_gff_path = None
libraries = None
sum_fc_thres = 100
base_path = None
variant = None

i = 1
while i < len(sys.argv):
    arg = sys.argv[i]
    if arg in ("--new-genome",):
        variant = "new_genome"
        i += 1
        continue
    if arg == "--variant":
        variant = sys.argv[i + 1]
        i += 2
        continue
    if arg == "--base-path":
        base_path = sys.argv[i + 1]
        i += 2
        continue
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
    elif arg == '--fc-pre-mirdeep':
        featurecounts_pre_mirdeep_path = sys.argv[i + 1]
    elif arg == '--fc-pre-sRNAbench':
        featurecounts_pre_sRNAbench_path = sys.argv[i + 1]
    elif arg in ('--fc-mirbase', '--fc_mirbase'):
        featurecounts_mirbase_path = sys.argv[i + 1]
    elif arg == '-rm':
        remaining_mirdeep_path = sys.argv[i + 1]
    elif arg == '-r1m':
        remaining1_mirdeep_path = sys.argv[i + 1]
    elif arg == '-r2m':
        remaining2_mirdeep_path = sys.argv[i + 1]
    elif arg == '-rs':
        remaining_sRNAbench_path = sys.argv[i + 1]
    elif arg == '-mgff':
        mirbase_gff_path = sys.argv[i + 1]
    elif arg == '-l':
        libraries = sys.argv[i + 1].split(',')
    elif arg == '--sum-fc-thres':
        sum_fc_thres = int(sys.argv[i + 1])
    elif arg in ('--mirdeep-mirbase-inter', '--mirdeep-mibrase-inter'):
        mirdeep_mibrase_inter = sys.argv[i + 1]
    elif arg == '--mirdeep-mirgenedb-inter':
        mirdeep_mirgenedb_inter = sys.argv[i + 1]
    elif arg in ('--sRNAbench-mirbase-inter', '--sRNAbench-mibrase-inter'):
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
              f' --new-genome: use new genome folder structure (Species_newGenome).\n'
              f' --mirdeep-inter-table <path>: path to bedtools -a mirdeep and -b sRNAbench intersection .bed file.\n'
              f' --sRNAbench-inter-table <path>: path to bedtools -a sRNAbench and -b mirdeep intersection .bed file.\n'
              f' --blast-mirdeep <path>: path to mirdeep blast results file.\n'
              f' --blast-sRNAbench <path>: path to sRNAbench blast results file.\n'
              f' --fc-mirdeep <path>: path to mirdeep featurecounts results file for mature/star (full counts, not the summary file).\n'
              f' --fc-sRNAbench <path>: path to sRNAbench featurecounts results file for mature/star (full counts, not the summary file).\n'
              f' --fc-pre-mirdeep <path>: path to mirdeep featurecounts results file for precursors (full counts, not the summary file).\n'
              f' --fc-pre-sRNAbench <path>: path to sRNAbench featurecounts results file for precursors (full counts, not the summary file).\n'
              f' -rm <path>: united mirdeep remaining file (mirdeep_all_remaining_filtered.csv).\n'
              f' -r1m/-r2m <path>: deprecated split remaining files (use -rm when possible).\n'
              f' -rs <path>: path to remaining sRNAbench candidates file, sRNAbench_remaining.csv.\n'
              f' -l <list>: list of sequencing libraries. Write the list seperated with commas, witout spaces. Example: library1,library2,library3 \n'
              f' --sum-fc-thres <int>: filtering threshold. Any candidates with sum_fc_m <= threshold will be filtered.\n'
              f'\nElegans only parameters:\n'
              f' --mirdeep-mirbase-inter <path>: path to bedtools -a mirdeep and -b mirbase intersection .bed file.\n'
              f' --mirdeep-mirgenedb-inter <path>: path to bedtools -a mirdeep and -b mirgenedb intersection .bed file.\n'
              f' --sRNAbench-mirbase-inter <path>: path to bedtools -a sRNAbench and -b mirbase intersection .bed file.\n'
              f' --sRNAbench-mirgenedb-inter <path>: path to bedtools -a sRNAbench and -b mirgenedb intersection .bed file.\n'
              f' --mirbase-mirgenedb-inter <path>: path to bedtools -a mirbase and -b mirgenedb intersection .bed file.\n'
              f' --mirbase-mirdeep-inter <path>: path to bedtools -a mirbase and -b mirdeep intersection .bed file.\n'
              f' --mirbase-sRNAbench-inter <path>: path to bedtools -a mirbase and -b sRNAbench intersection .bed file.\n'
              f' --fc-mirbase <path>: path to mirbase featurecounts results file (full counts, not the summary file).\n'
              )
        sys.exit()
    i += 2

cfg = get_species_config(species, base_path, variant=variant)
output_dir = cfg["output_dir"]
use_blast = species_uses_blast(cfg)
use_mirbase = species_uses_mirbase(cfg)

if use_mirbase:
    elegans_inputs = {
        "--mirdeep-mirbase-inter": mirdeep_mibrase_inter,
        "--mirdeep-mirgenedb-inter": mirdeep_mirgenedb_inter,
        "--sRNAbench-mirbase-inter": sRNAbench_mibrase_inter,
        "--sRNAbench-mirgenedb-inter": sRNAbench_mirgenedb_inter,
        "--mirbase-mirgenedb-inter": mirbase_mirgenedb_inter,
        "--mirbase-mirdeep-inter": mirbase_mirdeep_inter,
        "--mirbase-sRNAbench-inter": mirbase_sRNAbench_inter,
        "--fc-mirbase": featurecounts_mirbase_path,
        "-mgff": mirbase_gff_path,
    }
    missing = [flag for flag, value in elegans_inputs.items() if not value]
    if missing:
        print(
            "Error: {} requires these Elegans-specific arguments: {}".format(
                species, ", ".join(missing)
            ),
            file=sys.stderr,
        )
        sys.exit(2)

if remaining_mirdeep_path is None and remaining1_mirdeep_path and remaining2_mirdeep_path:
    r1 = pd.read_csv(remaining1_mirdeep_path, sep="\t")
    r2 = pd.read_csv(remaining2_mirdeep_path, sep="\t", names=r1.columns, skiprows=1)
    remaining_mirdeep = pd.concat([r1, r2], ignore_index=True)
else:
    remaining_mirdeep = pd.read_csv(remaining_mirdeep_path, sep="\t")

# -----CREATE INTERSECTIONS TABLES-----
# -----mirdeep intersections table:-----

mirdeep_intersection_columns = ['Chr_mirdeep', '.1', 'pre_miRNA1', 'Start_mirdeep', 'End_mirdeep', '.2', 'Strand_mirdeep', '.3', 'Description_mirdeep', 'Chr_sRNAbench', '.4', 'pre_miRNA2', 'Start_sRNAbench', 'End_sRNAbench', '.5', 'Strand_sRNAbench', '.6', 'Description_sRNAbench']
mirdeep_intersections_table = load_cross_intersections(
    mirdeep_intersections_table_path,
    mirdeep_intersection_columns,
    "mirdeep",
    "sRNAbench",
)
print(f"INITIAL SHAPE of mirdeep_intersections_table: {mirdeep_intersections_table.shape}")
print(f"Shape before deduplication: {mirdeep_intersections_table.shape}")

# --- FIX: Remove duplicate rows based on the unique description column ---
mirdeep_intersections_table.drop_duplicates(subset=['Description_mirdeep'], keep='first', inplace=True)

print(f"Shape after deduplication: {mirdeep_intersections_table.shape}")
mirdeep_intersections_table = mirdeep_intersections_table.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
n_srna_hit = int((mirdeep_intersections_table['Description_sRNAbench'] != '.').sum())
print(f"miRdeep rows with sRNAbench overlap: {n_srna_hit}/{len(mirdeep_intersections_table)}")
mirdeep_intersections_table['index'] = mirdeep_intersections_table['Description_mirdeep'].str.split(';').apply(lambda x: x[3]).str.replace('ID=','')


if use_mirbase:
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
        int)

    mirdeep_sRNAbench_mirbase = pd.merge(mirdeep_intersections_table, mirdeep_mirbase.iloc[:, 4:10], on='Description_mirdeep',
                                         how='left')
    mirdeep_sRNAbench_mirbase['T/F_mirbase'] = (mirdeep_sRNAbench_mirbase['Description_mirbase'] != '.').astype(
        int)
    mirdeep_intersections_table = pd.merge(mirdeep_sRNAbench_mirbase, mirdeep_mirgenedb.iloc[:, 4:10],
                                           on='Description_mirdeep', how='left')
    mirdeep_intersections_table['T/F_mirgenedb'] = (mirdeep_intersections_table['Description_mirgenedb'] != '.').astype(
        int)

# -----sRNAbench intersections table:-----

sRNAbench_intersection_columns = ['Chr_sRNAbench', '.1', 'pre_miRNA1', 'Start_sRNAbench', 'End_sRNAbench', '.2', 'Strand_sRNAbench', '.3', 'Description_sRNAbench', 'Chr_mirdeep', '.4', 'pre_miRNA2', 'Start_mirdeep', 'End_mirdeep', '.5', 'Strand_mirdeep', '.6', 'Description_mirdeep']
sRNAbench_intersections_table = load_cross_intersections(
    sRNAbench_intersections_table_path,
    sRNAbench_intersection_columns,
    "sRNAbench",
    "mirdeep",
)
sRNAbench_intersections_table = sRNAbench_intersections_table.drop(['.1', 'pre_miRNA1', '.2', '.3', '.4', 'pre_miRNA2', '.5', '.6'], axis=1)
sRNAbench_intersections_table['index'] = sRNAbench_intersections_table['Description_sRNAbench'].str.split(';').apply(lambda x: x[3]).str.replace('ID=','')


if use_mirbase:
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
        int)

    sRNAbench_mirdeep_mirbase = pd.merge(sRNAbench_intersections_table, sRNAbench_mirbase.iloc[:, 4:10], on='Description_sRNAbench',
                                         how='left')
    sRNAbench_mirdeep_mirbase['T/F_mirbase'] = (sRNAbench_mirdeep_mirbase['Description_mirbase'] != '.').astype(
        int)
    sRNAbench_intersections_table = pd.merge(sRNAbench_mirdeep_mirbase, sRNAbench_mirgenedb.iloc[:, 4:10],
                                             on='Description_sRNAbench', how='left')
    sRNAbench_intersections_table['T/F_mirgenedb'] = (
                sRNAbench_intersections_table['Description_mirgenedb'] != '.').astype(int)


if use_mirbase:
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
        int)

    mirbase_mirgenedb_mirdeep = pd.merge(mirbase_mirgenedb, mirbase_mirdeep.iloc[:, 4:10], on='Description_mirbase',
                                         how='left')
    mirbase_mirgenedb_mirdeep['T/F_mirdeep'] = (mirbase_mirgenedb_mirdeep['Description_mirdeep'] != '.').astype(
        int)
    mirbase_intersections_table = pd.merge(mirbase_mirgenedb_mirdeep, mirbase_sRNAbench.iloc[:, 4:10],
                                           on='Description_mirbase', how='left')
    mirbase_intersections_table['T/F_sRNAbench'] = (mirbase_intersections_table['Description_sRNAbench'] != '.').astype(
        int)

# -----ADD BLAST RESULTS (optional)-----
blast_mirdeep_orig = None
blast_sRNAbench_orig = None
if use_blast:
    if not blast_mirdeep_path or not blast_sRNAbench_path:
        print("Error: --blast-mirdeep and --blast-sRNAbench are required for this species.", file=sys.stderr)
        sys.exit(1)
    blast_mirdeep_orig, blast_mirdeep = load_blast_results(blast_mirdeep_path, "mirdeep")
    mirdeep_intersections_table = merge_blast(
        mirdeep_intersections_table, blast_mirdeep, "Description_mirdeep", index_pos=3
    )
    blast_sRNAbench_orig, blast_sRNAbench = load_blast_results(blast_sRNAbench_path, "srnabench")
    sRNAbench_intersections_table = merge_blast(
        sRNAbench_intersections_table, blast_sRNAbench, "Description_sRNAbench", index_pos=3
    )

# -----ADD FEATURE COUNTS:-----
# ---miRDeep:
featurecounts_mirdeep = pd.read_csv(featurecounts_mirdeep_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
featurecounts_mirdeep = featurecounts_mirdeep.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
featurecounts_mirdeep = featurecounts_mirdeep.iloc[2:]


# Create mature/star and index column
featurecounts_mirdeep['index'] = featurecounts_mirdeep['Geneid'].str.split('|').apply(lambda x : x[4])
featurecounts_mirdeep['mature/star'] = featurecounts_mirdeep['Geneid'].str.split('|').apply(lambda x : x[2])
featurecounts_mirdeep = featurecounts_mirdeep.drop('Geneid', axis=1)

# Casting libraries columns to int64
cast_dict = {k: 'int64' for k in libraries}
featurecounts_mirdeep = featurecounts_mirdeep.astype(cast_dict)

# Separate df into mature and star
mature_counts = featurecounts_mirdeep[featurecounts_mirdeep['mature/star'] == 'm'].copy()
libraries_mature = [library + '_m' for library in libraries]
rename_dict = dict(zip(libraries, libraries_mature))
mature_counts = mature_counts.rename(columns=rename_dict)
mature_counts['sum_FC_m'] = mature_counts[libraries_mature].sum(axis=1)
mature_counts = mature_counts.drop('mature/star', axis=1)


star_counts = featurecounts_mirdeep[featurecounts_mirdeep['mature/star'] == 's'].copy()
libraries_star = [library + '_s' for library in libraries]
rename_dict_star = dict(zip(libraries, libraries_star))
star_counts = star_counts.rename(columns=rename_dict_star)
star_counts['sum_FC_s'] = star_counts[libraries_star].sum(axis=1)
star_counts['sum_FC_s > 100?'] = np.where(star_counts['sum_FC_s'] > 100, 1, 0)
star_counts = star_counts.drop('mature/star', axis=1)

# Load precursor counts from its own file
featurecounts_pre_mirdeep = pd.read_csv(featurecounts_pre_mirdeep_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
featurecounts_pre_mirdeep = featurecounts_pre_mirdeep.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
featurecounts_pre_mirdeep = featurecounts_pre_mirdeep.iloc[2:]
featurecounts_pre_mirdeep["index"] = ["index={}".format(i) for i in range(len(featurecounts_pre_mirdeep))]
precursor_counts = featurecounts_pre_mirdeep


libraries_pre = [library + '_pre' for library in libraries]
rename_dict_pre = dict(zip(libraries, libraries_pre))
precursor_counts = precursor_counts.rename(columns=rename_dict_pre)
precursor_counts = precursor_counts.astype({k: 'int64' for k in libraries_pre})
precursor_counts['sum_FC_pre'] = precursor_counts[libraries_pre].sum(axis=1)


# Merge mirdeep & blast results and featurecounts results
mirdeep_blast_m_intersections_table = pd.merge(mirdeep_intersections_table, mature_counts, on='index', how='left')
mirdeep_blast_ms_intersections_table = pd.merge(mirdeep_blast_m_intersections_table, star_counts, on='index', how='left')
mirdeep_blast_fc_intersections_table = pd.merge(mirdeep_blast_ms_intersections_table, precursor_counts, on='index', how='left')
mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table.drop('index', axis=1)

# filter by sum_fc_m < threshold
mirdeep_blast_fc_intersections_table["sum_FC_m > thres"] = np.where(mirdeep_blast_fc_intersections_table['sum_FC_m'] > sum_fc_thres, 1, 0)

# Extract readcounts columns (partner Description may be '.' / NaN under bedtools -loj)
mirdeep_blast_fc_intersections_table['RC_m mirdeep'] = extract_description_field(mirdeep_blast_fc_intersections_table["Description_mirdeep"], 1)
mirdeep_blast_fc_intersections_table['RC_s mirdeep'] = extract_description_field(mirdeep_blast_fc_intersections_table["Description_mirdeep"], 2)
mirdeep_blast_fc_intersections_table['RC_m sRNAbench'] = extract_description_field(mirdeep_blast_fc_intersections_table["Description_sRNAbench"], 1)
mirdeep_blast_fc_intersections_table['RC_s sRNAbench'] = extract_description_field(mirdeep_blast_fc_intersections_table["Description_sRNAbench"], 2)
mirdeep_blast_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']] = mirdeep_blast_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']].fillna('0')
mirdeep_blast_fc_intersections_table['RC_m mirdeep'] = mirdeep_blast_fc_intersections_table['RC_m mirdeep'].astype(str).str.replace('RC_m=', '').astype('int64')
mirdeep_blast_fc_intersections_table['RC_s mirdeep'] = mirdeep_blast_fc_intersections_table['RC_s mirdeep'].astype(str).str.replace('RC_s=', '').astype('int64')
mirdeep_blast_fc_intersections_table['RC_m sRNAbench'] = mirdeep_blast_fc_intersections_table['RC_m sRNAbench'].astype(str).str.replace('RC_m=', '').astype('int64')
mirdeep_blast_fc_intersections_table['RC_s sRNAbench'] = mirdeep_blast_fc_intersections_table['RC_s sRNAbench'].astype(str).str.replace('RC_s=', '').astype('int64')

# Create diff columns
mirdeep_blast_fc_intersections_table['Diff Sum_FC_m / RC_m mirdeep'] = mirdeep_blast_fc_intersections_table['sum_FC_m'] / mirdeep_blast_fc_intersections_table['RC_m mirdeep']
mirdeep_blast_fc_intersections_table['Diff Sum_FC_m / RC_m sRNAbench'] = mirdeep_blast_fc_intersections_table['sum_FC_m'] / mirdeep_blast_fc_intersections_table['RC_m sRNAbench']
mirdeep_blast_fc_intersections_table['Diff Sum_FC_s / RC_s mirdeep'] = mirdeep_blast_fc_intersections_table['sum_FC_s'] / mirdeep_blast_fc_intersections_table['RC_s mirdeep']
mirdeep_blast_fc_intersections_table['Diff Sum_FC_s / RC_s sRNAbench'] = mirdeep_blast_fc_intersections_table['sum_FC_s'] / mirdeep_blast_fc_intersections_table['RC_s sRNAbench']

# Normalize featurecounts to reads per million
mature_rpm = [column + "_rpm" for column in libraries_mature]
star_rpm = [column + "_rpm" for column in libraries_star]
pre_rpm = [column + "_rpm" for column in libraries_pre]

for i in range(0, len(libraries)):
    library_m = libraries_mature[i]
    library_s = libraries_star[i]
    library_pre = libraries_pre[i]
    total = mirdeep_blast_fc_intersections_table[[library_m, library_s, library_pre]].sum().sum()
    mirdeep_blast_fc_intersections_table[[mature_rpm[i], star_rpm[i], pre_rpm[i]]] = round((mirdeep_blast_fc_intersections_table[[library_m, library_s, library_pre]] / total) * 1000000, 0)


mirdeep_blast_fc_intersections_table['sum_FC_m_rpm'] = mirdeep_blast_fc_intersections_table[mature_rpm].sum(axis=1)
mirdeep_blast_fc_intersections_table['sum_FC_s_rpm'] = mirdeep_blast_fc_intersections_table[star_rpm].sum(axis=1)
mirdeep_blast_fc_intersections_table['sum_FC_pre_rpm'] = mirdeep_blast_fc_intersections_table[pre_rpm].sum(axis=1)
mirdeep_blast_fc_intersections_table['mean_m_rpm'] = round(mirdeep_blast_fc_intersections_table[mature_rpm].mean(axis=1), 2)
mirdeep_blast_fc_intersections_table['mean_s_rpm'] = round(mirdeep_blast_fc_intersections_table[star_rpm].mean(axis=1), 2)
mirdeep_blast_fc_intersections_table['mean_pre_rpm'] = round(mirdeep_blast_fc_intersections_table[pre_rpm].mean(axis=1), 2)

# MODIFIED: Calculate mature/precursor read count ratios for mirdeep
ratio_columns = []
for library in libraries:
    mature_col = f'{library}_m'
    pre_col = f'{library}_pre'
    ratio_col = f'{library}_m/pre_ratio'
    ratio_columns.append(ratio_col)
    # Use numpy to handle division by zero safely
    mirdeep_blast_fc_intersections_table[ratio_col] = np.divide(
        mirdeep_blast_fc_intersections_table[mature_col],
        mirdeep_blast_fc_intersections_table[pre_col]
    ).fillna(0)
mirdeep_blast_fc_intersections_table.replace([np.inf, -np.inf], 0, inplace=True)

# ---sRNAbench:
featurecounts_sRNAbench = pd.read_csv(featurecounts_sRNAbench_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
featurecounts_sRNAbench = featurecounts_sRNAbench.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
featurecounts_sRNAbench = featurecounts_sRNAbench.iloc[2:]

# Create 5p/3p column
featurecounts_sRNAbench['5p/3p'] = featurecounts_sRNAbench['Geneid'].str.split('|', expand=True)[0]
featurecounts_sRNAbench['5p/3p'] = featurecounts_sRNAbench['5p/3p'].str.split('-').apply(lambda x : x[-1])
featurecounts_sRNAbench['5p/3p'] = featurecounts_sRNAbench['5p/3p'].str.split('_', expand=True)[0]
featurecounts_sRNAbench = featurecounts_sRNAbench.rename(columns={'5p/3p':'mature'})

# Create mature/star and index column
featurecounts_sRNAbench['index'] = featurecounts_sRNAbench['Geneid'].str.split('|').apply(lambda x : x[3])
featurecounts_sRNAbench['mature/star'] = featurecounts_sRNAbench['Geneid'].str.split('|').apply(lambda x : x[1])
featurecounts_sRNAbench = featurecounts_sRNAbench.drop('Geneid', axis=1)

# Casting libraries columns to int64
featurecounts_sRNAbench = featurecounts_sRNAbench.astype(cast_dict)

# Separate df into mature and star
mature_counts = featurecounts_sRNAbench[featurecounts_sRNAbench['mature/star'] == 'm'].copy()
mature_counts = mature_counts.rename(columns=rename_dict)
mature_counts['sum_FC_m'] = mature_counts[libraries_mature].sum(axis=1)
mature_counts = mature_counts.drop('mature/star', axis=1)

star_counts = featurecounts_sRNAbench[featurecounts_sRNAbench['mature/star'] == 's'].copy()
star_counts = star_counts.rename(columns=rename_dict_star)
star_counts['sum_FC_s'] = star_counts[libraries_star].sum(axis=1)
star_counts['sum_FC_s > 100?'] = np.where(star_counts['sum_FC_s'] > 100, 1, 0)
star_counts = star_counts.drop(['mature/star', 'mature'], axis=1)

# Load precursor counts from its own file for sRNAbench
featurecounts_pre_sRNAbench = pd.read_csv(featurecounts_pre_sRNAbench_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
featurecounts_pre_sRNAbench = featurecounts_pre_sRNAbench.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
featurecounts_pre_sRNAbench = featurecounts_pre_sRNAbench.iloc[2:]
featurecounts_pre_sRNAbench["index"] = ["index={}".format(i) for i in range(len(featurecounts_pre_sRNAbench))]
precursor_counts_srna = featurecounts_pre_sRNAbench
precursor_counts_srna = precursor_counts_srna.rename(columns=rename_dict_pre)
precursor_counts_srna = precursor_counts_srna.astype({k: 'int64' for k in libraries_pre})
precursor_counts_srna['sum_FC_pre'] = precursor_counts_srna[libraries_pre].sum(axis=1)


# Merge sRNAbench & blast results and featurecounts results
sRNAbench_blast_m_intersections_table = pd.merge(sRNAbench_intersections_table, mature_counts, on='index', how='left')
sRNAbench_blast_ms_intersections_table = pd.merge(sRNAbench_blast_m_intersections_table, star_counts, on='index', how='left')
sRNAbench_blast_fc_intersections_table = pd.merge(sRNAbench_blast_ms_intersections_table, precursor_counts_srna, on='index', how='left')
sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table.drop('index', axis=1)


# filter by sum_fc_m < threshold
sRNAbench_blast_fc_intersections_table["sum_FC_m > thres"] = np.where(sRNAbench_blast_fc_intersections_table['sum_FC_m'] > sum_fc_thres, 1, 0)

# Extract readcounts columns (partner Description may be '.' / NaN under bedtools -loj)
sRNAbench_blast_fc_intersections_table['RC_m mirdeep'] = extract_description_field(sRNAbench_blast_fc_intersections_table["Description_mirdeep"], 1)
sRNAbench_blast_fc_intersections_table['RC_s mirdeep'] = extract_description_field(sRNAbench_blast_fc_intersections_table["Description_mirdeep"], 2)
sRNAbench_blast_fc_intersections_table['RC_m sRNAbench'] = extract_description_field(sRNAbench_blast_fc_intersections_table["Description_sRNAbench"], 1)
sRNAbench_blast_fc_intersections_table['RC_s sRNAbench'] = extract_description_field(sRNAbench_blast_fc_intersections_table["Description_sRNAbench"], 2)
sRNAbench_blast_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']] = sRNAbench_blast_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']].fillna('0')
sRNAbench_blast_fc_intersections_table['RC_m mirdeep'] = sRNAbench_blast_fc_intersections_table['RC_m mirdeep'].astype(str).str.replace('RC_m=', '').astype('int64')
sRNAbench_blast_fc_intersections_table['RC_s mirdeep'] = sRNAbench_blast_fc_intersections_table['RC_s mirdeep'].astype(str).str.replace('RC_s=', '').astype('int64')
sRNAbench_blast_fc_intersections_table['RC_m sRNAbench'] = sRNAbench_blast_fc_intersections_table['RC_m sRNAbench'].astype(str).str.replace('RC_m=', '').astype('int64')
sRNAbench_blast_fc_intersections_table['RC_s sRNAbench'] = sRNAbench_blast_fc_intersections_table['RC_s sRNAbench'].astype(str).str.replace('RC_s=', '').astype('int64')

# Create diff columns
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_m / RC_m mirdeep'] = sRNAbench_blast_fc_intersections_table['sum_FC_m'] / sRNAbench_blast_fc_intersections_table['RC_m mirdeep']
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_m / RC_m sRNAbench'] = sRNAbench_blast_fc_intersections_table['sum_FC_m'] / sRNAbench_blast_fc_intersections_table['RC_m sRNAbench']
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_s / RC_s mirdeep'] = sRNAbench_blast_fc_intersections_table['sum_FC_s'] / sRNAbench_blast_fc_intersections_table['RC_s mirdeep']
sRNAbench_blast_fc_intersections_table['Diff Sum_FC_s / RC_s sRNAbench'] = sRNAbench_blast_fc_intersections_table['sum_FC_s'] / sRNAbench_blast_fc_intersections_table['RC_s sRNAbench']

# Normalize featurecounts to reads per million
for i in range(0, len(libraries)):
    library_m = libraries_mature[i]
    library_s = libraries_star[i]
    library_pre = libraries_pre[i]
    total = sRNAbench_blast_fc_intersections_table[[library_m, library_s, library_pre]].sum().sum()
    sRNAbench_blast_fc_intersections_table[[mature_rpm[i], star_rpm[i], pre_rpm[i]]] = round((sRNAbench_blast_fc_intersections_table[[library_m, library_s, library_pre]] / total) * 1000000, 0)


sRNAbench_blast_fc_intersections_table['sum_FC_m_rpm'] = sRNAbench_blast_fc_intersections_table[mature_rpm].sum(axis=1)
sRNAbench_blast_fc_intersections_table['sum_FC_s_rpm'] = sRNAbench_blast_fc_intersections_table[star_rpm].sum(axis=1)
sRNAbench_blast_fc_intersections_table['sum_FC_pre_rpm'] = sRNAbench_blast_fc_intersections_table[pre_rpm].sum(axis=1)
sRNAbench_blast_fc_intersections_table['mean_m_rpm'] = round(sRNAbench_blast_fc_intersections_table[mature_rpm].mean(axis=1), 2)
sRNAbench_blast_fc_intersections_table['mean_s_rpm'] = round(sRNAbench_blast_fc_intersections_table[star_rpm].mean(axis=1), 2)
sRNAbench_blast_fc_intersections_table['mean_pre_rpm'] = round(sRNAbench_blast_fc_intersections_table[pre_rpm].mean(axis=1), 2)

# MODIFIED: Calculate mature/precursor read count ratios for sRNAbench
ratio_columns_srna = []
for library in libraries:
    mature_col = f'{library}_m'
    pre_col = f'{library}_pre'
    ratio_col = f'{library}_m/pre_ratio'
    ratio_columns_srna.append(ratio_col)
    # Use numpy to handle division by zero safely
    sRNAbench_blast_fc_intersections_table[ratio_col] = np.divide(
        sRNAbench_blast_fc_intersections_table[mature_col],
        sRNAbench_blast_fc_intersections_table[pre_col]
    ).fillna(0)
sRNAbench_blast_fc_intersections_table.replace([np.inf, -np.inf], 0, inplace=True)


if use_mirbase:
    # -----miRBase:
    featurecounts_mirbase = pd.read_csv(featurecounts_mirbase_path, sep='\t', names=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'] + libraries)
    featurecounts_mirbase = featurecounts_mirbase.drop(['Chr', 'Start', 'End', 'Strand', 'Length'], axis=1)
    featurecounts_mirbase = featurecounts_mirbase.iloc[2:]


    # Create index column for featurecounts
    featurecounts_mirbase['index'] = featurecounts_mirbase['Geneid'].str.split(';').apply(lambda x : x[3])
    featurecounts_mirbase['index'] = featurecounts_mirbase['index'].str.replace('Derives_from=MI', '')

    # Create 5p/3p columns
    featurecounts_mirbase['5p/3p'] = featurecounts_mirbase['Geneid'].str.split(';').apply(lambda x : x[2])
    featurecounts_mirbase['5p/3p'] = featurecounts_mirbase['5p/3p'].str[-2:]

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

    # Create empty mature df and star df, iterate the rows of 5p and 3p and add to the relevant.
    mature_df = pd.DataFrame(columns=libraries_mature + ['mature'])
    star_df = pd.DataFrame(columns=libraries_star)
    for i in range(0, len(counts_mb_5p)):
        row_5p = counts_mb_5p.iloc[i]
        row_3p = counts_mb_3p.iloc[i]
        if row_5p['sum'] > row_3p['sum']: # if 5p is the mature
            rename_dict = dict(zip(libraries_5p, libraries_mature))
            row_5p = row_5p.rename(rename_dict)
            rename_dict = dict(zip(libraries_3p, libraries_star))
            row_3p = row_3p.rename(rename_dict)
            row_5p['mature'] = '5p'
            mature_df = mature_df.append(row_5p)
            star_df = star_df.append(row_3p)
        elif row_5p['sum'] <= row_3p['sum']: # else 3p is the mature
            rename_dict = dict(zip(libraries_5p, libraries_star))
            row_5p = row_5p.rename(rename_dict)
            rename_dict = dict(zip(libraries_3p, libraries_mature))
            row_3p = row_3p.rename(rename_dict)
            row_3p['mature'] = '3p'
            mature_df = mature_df.append(row_3p)
            star_df = star_df.append(row_5p)

    mature_df = mature_df.rename(columns={'sum':'sum_FC_m'})
    star_df = star_df.rename(columns={'sum':'sum_FC_s'})
    star_df['sum_FC_s > 100?'] = np.where(star_df['sum_FC_s'] > 100, 1, 0)

    # Those that are only one strand and are no 5p/3p are determined as mature.
    counts_no_5p3p = counts_no_5p3p.drop('5p/3p', axis=1)
    rename_dict = dict(zip(libraries, libraries_mature))
    counts_no_5p3p = counts_no_5p3p.rename(columns=rename_dict)
    counts_no_5p3p['mature'] = counts_no_5p3p['Geneid'].str.split(';', expand=True)[5]
    counts_no_5p3p['sum_FC_m'] = counts_no_5p3p[libraries_mature].sum(axis=1)
    mature_df = mature_df.append(counts_no_5p3p)

    counts_no_5p3p_filler = counts_no_5p3p.copy()
    rename_dict = dict(zip(libraries_mature, libraries_star))
    counts_no_5p3p_filler = counts_no_5p3p_filler.rename(columns=rename_dict)
    counts_no_5p3p_filler = counts_no_5p3p_filler.drop(['sum_FC_m', 'mature'], axis=1)
    counts_no_5p3p_filler[libraries_star] = 0
    counts_no_5p3p_filler['sum_FC_s'] = 0
    counts_no_5p3p_filler['sum_FC_s > 100?'] = 0
    star_df = star_df.append(counts_no_5p3p_filler)

    # Create index column for mirbase
    mirbase_intersections_table['index'] = mirbase_intersections_table['Description_mirbase'].str.split(';').apply(lambda x : x[0])
    mirbase_intersections_table['index'] = mirbase_intersections_table['index'].str.replace('ID=MI', '')

    # Merge mirbase results and mirbase featurecounts results
    mirbase_m_intersections_table = pd.merge(mirbase_intersections_table, mature_df, on='index', how='left')
    mirbase_fc_intersections_table = pd.merge(mirbase_m_intersections_table, star_df, on='index', how='left')
    mirbase_fc_intersections_table = mirbase_fc_intersections_table.drop('index', axis=1)

    # filter by sum_fc_m < threshold
    mirbase_fc_intersections_table["sum_FC_m > thres"] = np.where(mirbase_fc_intersections_table['sum_FC_m'] > sum_fc_thres, 1, 0)

    # Extract readcounts columns (partner Description may be '.' / NaN under bedtools -loj)
    mirbase_fc_intersections_table['RC_m mirdeep'] = extract_description_field(mirbase_fc_intersections_table["Description_mirdeep"], 1)
    mirbase_fc_intersections_table['RC_s mirdeep'] = extract_description_field(mirbase_fc_intersections_table["Description_mirdeep"], 2)
    mirbase_fc_intersections_table['RC_m sRNAbench'] = extract_description_field(mirbase_fc_intersections_table["Description_sRNAbench"], 1)
    mirbase_fc_intersections_table['RC_s sRNAbench'] = extract_description_field(mirbase_fc_intersections_table["Description_sRNAbench"], 2)
    mirbase_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']] = mirbase_fc_intersections_table[['RC_m mirdeep', 'RC_s mirdeep']].fillna('0')
    mirbase_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']] = mirbase_fc_intersections_table[['RC_m sRNAbench', 'RC_s sRNAbench']].fillna('0')
    mirbase_fc_intersections_table['RC_m mirdeep'] = mirbase_fc_intersections_table['RC_m mirdeep'].astype(str).str.replace('RC_m=', '').astype('int64')
    mirbase_fc_intersections_table['RC_s mirdeep'] = mirbase_fc_intersections_table['RC_s mirdeep'].astype(str).str.replace('RC_s=', '').astype('int64')
    mirbase_fc_intersections_table['RC_m sRNAbench'] = mirbase_fc_intersections_table['RC_m sRNAbench'].astype(str).str.replace('RC_m=', '').astype('int64')
    mirbase_fc_intersections_table['RC_s sRNAbench'] = mirbase_fc_intersections_table['RC_s sRNAbench'].astype(str).str.replace('RC_s=', '').astype('int64')

    # Create diff columns
    mirbase_fc_intersections_table['Diff Sum_FC_m / RC_m mirdeep'] = mirbase_fc_intersections_table['sum_FC_m'] / mirbase_fc_intersections_table['RC_m mirdeep']
    mirbase_fc_intersections_table['Diff Sum_FC_m / RC_m sRNAbench'] = mirbase_fc_intersections_table['sum_FC_m'] / mirbase_fc_intersections_table['RC_m sRNAbench']
    mirbase_fc_intersections_table['Diff Sum_FC_s / RC_s mirdeep'] = mirbase_fc_intersections_table['sum_FC_s'] / mirbase_fc_intersections_table['RC_s mirdeep']
    mirbase_fc_intersections_table['Diff Sum_FC_s / RC_s sRNAbench'] = mirbase_fc_intersections_table['sum_FC_s'] / mirbase_fc_intersections_table['RC_s sRNAbench']

    # Normalize featurecounts to reads per million
    cast_dict_rpm = {k: 'int64' for k in libraries_mature + libraries_star}
    mirbase_fc_intersections_table = mirbase_fc_intersections_table.astype(cast_dict_rpm)
    for i in range(0, len(libraries)):
        library_m = libraries_mature[i]
        library_s = libraries_star[i]
        total = mirbase_fc_intersections_table[[library_m, library_s]].sum().sum()
        mirbase_fc_intersections_table[[mature_rpm[i], star_rpm[i]]] = round((mirbase_fc_intersections_table[[library_m, library_s]] / total) * 1000000, 0)

    mirbase_fc_intersections_table['sum_FC_m_rpm'] = mirbase_fc_intersections_table[mature_rpm].sum(axis=1)
    mirbase_fc_intersections_table['sum_FC_s_rpm'] = mirbase_fc_intersections_table[star_rpm].sum(axis=1)
    mirbase_fc_intersections_table['mean_m_rpm'] = round(mirbase_fc_intersections_table[mature_rpm].mean(axis=1), 2)
    mirbase_fc_intersections_table['mean_s_rpm'] = round(mirbase_fc_intersections_table[star_rpm].mean(axis=1), 2)

# -----ADD SEQUENCES-----
# ---miRdeep:
print(f"INITIAL SHAPE of remaining_mirdeep: {remaining_mirdeep.shape}")
mirdeep_blast_fc_intersections_table = attach_remaining_columns(
    mirdeep_blast_fc_intersections_table,
    remaining_mirdeep,
    "Description_mirdeep",
    {
        "consensus mature sequence": "consensus mature sequence",
        "consensus star sequence": "consensus star sequence",
        "consensus precursor sequence": "consensus precursor sequence",
        "overlaps": "overlaps",
    },
    "miRdeep",
)
for column in [
    "consensus mature sequence",
    "consensus star sequence",
    "consensus precursor sequence",
]:
    mirdeep_blast_fc_intersections_table[column] = (
        mirdeep_blast_fc_intersections_table[column].str.upper()
    )
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
sRNAbench_blast_fc_intersections_table = attach_remaining_columns(
    sRNAbench_blast_fc_intersections_table,
    remaining_sRNAbench,
    "Description_sRNAbench",
    {
        "5pseq": "5pseq",
        "3pseq": "3pseq",
        "hairpinSeq": "hairpinSeq",
        "overlaps": "overlaps",
    },
    "sRNAbench",
)
for column in ["5pseq", "3pseq", "hairpinSeq"]:
    sRNAbench_blast_fc_intersections_table[column] = (
        sRNAbench_blast_fc_intersections_table[column].str.replace('T', 'U')
    )

# Extract loop size
def loop_size(row):
    if (row["hairpinSeq"].find(str(row['5pseq'])) == -1) or (row["hairpinSeq"].find(str(row['3pseq'])) == -1):
        return -1
    index_end_5p = len((str(row['5pseq'])))
    index_start_3p = row["hairpinSeq"].find(str(row['3pseq']))
    loop_size = index_start_3p - index_end_5p
    return loop_size

sRNAbench_blast_fc_intersections_table['loop_size'] = sRNAbench_blast_fc_intersections_table.apply(lambda row : loop_size(row), axis=1)
sRNAbench_blast_fc_intersections_table['mature_size'] = np.where(sRNAbench_blast_fc_intersections_table['mature'] == '5p', sRNAbench_blast_fc_intersections_table['5pseq'].str.len(), sRNAbench_blast_fc_intersections_table['3pseq'].str.len())
sRNAbench_blast_fc_intersections_table['star_size'] = np.where(sRNAbench_blast_fc_intersections_table['mature'] == '5p', sRNAbench_blast_fc_intersections_table['3pseq'].str.len(), sRNAbench_blast_fc_intersections_table['5pseq'].str.len())

# ---Mirbase
if use_mirbase:
    mirbase_fc_intersections_table['hairpinSeq'] = mirbase_fc_intersections_table['Description_mirbase'].str.split(';', expand=True)[3]
    mirbase_fc_intersections_table['Derives_from'] = mirbase_fc_intersections_table['Description_mirbase'].str.split(';', expand=True)[0]
    mirbase_fc_intersections_table['Derives_from'] = mirbase_fc_intersections_table['Derives_from'].str.split('=', expand=True)[1]
    annotation = gffpd.read_gff3(mirbase_gff_path)
    gff = annotation.df
    gff = gff[gff['type'] == 'miRNA'].copy()
    gff['Derives_from'] = gff['attributes'].str.split(';', expand=True)[3]
    gff['Derives_from'] = gff['Derives_from'].str.split('=', expand=True)[1]
    gff['sequence'] = gff['attributes'].str.split(';', expand=True)[4]
    gff['5p/3p'] = gff['attributes'].str.split(';', expand=True)[5]
    gff.loc[gff['5p/3p'] == '5p', '5pseq'] = gff['sequence']
    gff.loc[gff['5p/3p'] == '3p', '3pseq'] = gff['sequence']
    seq5p = gff[['Derives_from', '5pseq']]
    seq3p = gff[['Derives_from', '3pseq']]
    seq5p = seq5p.dropna()
    seq3p = seq3p.dropna()
    seq5p3p = pd.merge(seq5p, seq3p, left_on='Derives_from', right_on='Derives_from', how='outer')
    mirbase_fc_intersections_table = pd.merge(mirbase_fc_intersections_table, seq5p3p, left_on='Derives_from', right_on='Derives_from', how='left')

    mirbase_fc_intersections_table['loop_size'] = mirbase_fc_intersections_table.apply(lambda row : loop_size(row), axis=1)
    mirbase_fc_intersections_table['mature_size'] = np.where(mirbase_fc_intersections_table['mature'] == '5p', mirbase_fc_intersections_table['5pseq'].str.len(), mirbase_fc_intersections_table['3pseq'].str.len())
    mirbase_fc_intersections_table['star_size'] = np.where(mirbase_fc_intersections_table['mature'] == '5p', mirbase_fc_intersections_table['3pseq'].str.len(), mirbase_fc_intersections_table['5pseq'].str.len())
    mirbase_fc_intersections_table['star_size'] = mirbase_fc_intersections_table['star_size'].fillna(0)

# -----REORDER COLUMNS:-----

if use_mirbase:
    elegans_columns_mirdeep = ['T/F_sRNAbench', 'Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase', 'T/F_mirbase',
                       'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb']
    elegans_columns_sRNAbench = ['T/F_mirdeep', 'Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase', 'T/F_mirbase',
                       'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb']
    mirbase_fc_intersections_table = mirbase_fc_intersections_table[['Chr_mirbase', 'Start_mirbase', 'End_mirbase', 'Strand_mirbase', 'Description_mirbase',
                                                                     'Chr_mirgenedb', 'Start_mirgenedb', 'End_mirgenedb', 'Strand_mirgenedb', 'Description_mirgenedb', 'T/F_mirgenedb',
                                                                     'Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep', 'T/F_mirdeep',
                                                                     'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench', 'T/F_sRNAbench'] +
                                                                    libraries_mature + ['sum_FC_m', 'sum_FC_m > thres', 'RC_m mirdeep', 'RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench'] +
                                                                    libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s mirdeep', 'RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench'] +
                                                                    mature_rpm + ['sum_FC_m_rpm', 'mean_m_rpm'] + star_rpm + ['sum_FC_s_rpm', 'mean_s_rpm'] +
                                                                    ['5pseq', '3pseq', 'hairpinSeq', 'mature', 'mature_size', 'star_size', 'loop_size']
                                                                    ]
else:
    elegans_columns_mirdeep = []
    elegans_columns_sRNAbench = []


mirdeep_blast_fc_intersections_table = mirdeep_blast_fc_intersections_table[['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep',
                                                                             'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench'] + elegans_columns_mirdeep +
                                                                             libraries_mature + ['sum_FC_m', 'sum_FC_m > thres', 'RC_m mirdeep', 'RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench'] +
                                                                             libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s mirdeep', 'RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench'] +
                                                                             libraries_pre + ['sum_FC_pre'] +
                                                                             ratio_columns +
                                                                             mature_rpm + ['sum_FC_m_rpm', 'mean_m_rpm'] + star_rpm + ['sum_FC_s_rpm', 'mean_s_rpm'] + pre_rpm + ['sum_FC_pre_rpm', 'mean_pre_rpm'] +
                                                                             ['consensus mature sequence', 'consensus star sequence', 'consensus precursor sequence', 'mature', 'mature_size', 'star_size', 'loop_size', 'overlaps']
                                                                            ]

sRNAbench_blast_fc_intersections_table = sRNAbench_blast_fc_intersections_table[['Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'Description_sRNAbench',
                                                                             'Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Description_mirdeep'] + elegans_columns_sRNAbench +
                                                                             libraries_mature + ['sum_FC_m', 'sum_FC_m > thres', 'RC_m sRNAbench', 'RC_m mirdeep', 'Diff Sum_FC_m / RC_m sRNAbench', 'Diff Sum_FC_m / RC_m mirdeep'] +
                                                                             libraries_star + ['sum_FC_s', 'sum_FC_s > 100?', 'RC_s sRNAbench', 'RC_s mirdeep', 'Diff Sum_FC_s / RC_s sRNAbench', 'Diff Sum_FC_s / RC_s mirdeep'] +
                                                                             libraries_pre + ['sum_FC_pre'] +
                                                                             ratio_columns_srna +
                                                                             mature_rpm + ['sum_FC_m_rpm', 'mean_m_rpm'] + star_rpm + ['sum_FC_s_rpm', 'mean_s_rpm'] + pre_rpm + ['sum_FC_pre_rpm', 'mean_pre_rpm'] +
                                                                             ['5pseq', '3pseq', 'hairpinSeq', 'mature', 'mature_size', 'star_size', 'loop_size', 'overlaps']
                                                                                ]

# ------ADD TYPES------

if use_mirbase:
    # ------miRdeep:
    mirdeep_blast_fc_intersections_table['Type'] = np.zeros(len(mirdeep_blast_fc_intersections_table))
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 1
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 2
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 5
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 6
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 1) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 9
    mirdeep_blast_fc_intersections_table.loc[((mirdeep_blast_fc_intersections_table['T/F_sRNAbench'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirbase'] == 0) & (mirdeep_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 10


    # ------sRNAbench:
    sRNAbench_blast_fc_intersections_table['Type'] = np.zeros(len(sRNAbench_blast_fc_intersections_table))
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 1
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 1)), 'Type'] = 3
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 5
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 7
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 1) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 9
    sRNAbench_blast_fc_intersections_table.loc[((sRNAbench_blast_fc_intersections_table['T/F_mirdeep'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirbase'] == 0) & (sRNAbench_blast_fc_intersections_table['T/F_mirgenedb'] == 0)), 'Type'] = 11

    # ------mirbase:
    mirbase_fc_intersections_table['Type'] = np.zeros(len(mirbase_fc_intersections_table))
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 1
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 5
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 2
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 3
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 1) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 6
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 1)), 'Type'] = 7
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 1) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 4
    mirbase_fc_intersections_table.loc[((mirbase_fc_intersections_table['T/F_mirgenedb'] == 0) & (mirbase_fc_intersections_table['T/F_mirdeep'] == 0) & (mirbase_fc_intersections_table['T/F_sRNAbench'] == 0)), 'Type'] = 8

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

unified = unified.drop(['Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench', 'consensus mature sequence', 'consensus star sequence'], axis=1)
if use_mirbase:
    only_sRNAbench = sRNAbench_blast_fc_intersections_table[sRNAbench_blast_fc_intersections_table['Type'].isin([3, 7, 11])]
    only_mirbase = mirbase_fc_intersections_table[mirbase_fc_intersections_table['Type'].isin([4, 8])]
else:
    only_sRNAbench = sRNAbench_blast_fc_intersections_table[sRNAbench_blast_fc_intersections_table['Type'] == 3]


only_sRNAbench = only_sRNAbench.drop(['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep'], axis=1)
unified = unified.rename(columns={'Chr_mirdeep':'Chr', 'Start_mirdeep':'Start', 'End_mirdeep':'End', 'Strand_mirdeep':'Strand', 'consensus precursor sequence':'hairpinSeq'})
only_sRNAbench = only_sRNAbench.rename(columns={'Chr_sRNAbench':'Chr', 'Start_sRNAbench':'Start', 'End_sRNAbench':'End', 'Strand_sRNAbench':'Strand'})
if use_mirbase:
    only_mirbase = only_mirbase.drop(['Chr_mirdeep', 'Start_mirdeep', 'End_mirdeep', 'Strand_mirdeep', 'Chr_sRNAbench', 'Start_sRNAbench', 'End_sRNAbench', 'Strand_sRNAbench'], axis=1)
    only_mirbase = only_mirbase.rename(columns={'Chr_mirbase':'Chr', 'Start_mirbase':'Start', 'End_mirbase':'End', 'Strand_mirbase':'Strand'})
    unified = pd.concat([unified, only_sRNAbench, only_mirbase], axis=0, ignore_index=True)
    unified = unified.drop(['T/F_mirdeep', 'T/F_sRNAbench', 'T/F_mirbase', 'T/F_mirgenedb'], axis=1)
else:
    unified = pd.concat([unified, only_sRNAbench], axis=0, ignore_index=True)
# reorder sequences columns
columns = list(unified.columns)
i1, i2, i3 = columns.index('mature'), columns.index('5pseq'), columns.index('3pseq')
columns[i1], columns[i2], columns[i3] = columns[i2], columns[i3], columns[i1]
# reorder sRNAbench and mirbase description columns
columns.remove('Description_sRNAbench')
columns.insert(columns.index('Description_mirdeep') + 1, 'Description_sRNAbench')
if use_mirbase:
    columns.remove('Description_mirbase')
    columns.insert(columns.index('Description_mirdeep') + 2, 'Description_mirbase')
unified = unified[columns]


# --- Extract seed
unified['Seed'] = np.where(unified["mature"] == '5p', unified["5pseq"].str[1:8], unified["3pseq"].str[1:8])

seed_families = load_seed_table(cfg["seed_path"], cfg["seed_encoding"], cfg.get("seed_sep", "\t"))
seed_families = seed_families[['Family', 'Seed']]
seed_families = seed_families.drop_duplicates(subset='Seed')
unified = pd.merge(unified, seed_families, left_on='Seed', right_on='Seed', how='left')
unified['Family'].fillna(' ', inplace=True)
unified['Family'] = unified['Family'].str.replace(' ', 'UNKNOWN')

unified = unified.dropna(axis=1, thresh=1)
unified = unified.reindex(columns=[col for col in unified.columns if col != 'Type'] + ['Type'])

# --- Removing novel451
unified["novel451"] = np.where(
    normalize_description(unified['Description_sRNAbench']).str.contains("novel451"), 1, 0
)

# -----SAVE TO EXCEL-----

writer = pd.ExcelWriter(output_dir + 'intersections_table_{}.xlsx'.format(species))
mirdeep_blast_fc_intersections_table.to_excel(writer, sheet_name='miRdeep')
sRNAbench_blast_fc_intersections_table.to_excel(writer, sheet_name='sRNAbench')
if use_mirbase:
    mirbase_fc_intersections_table.to_excel(writer, sheet_name='mirbase')
unified.to_excel(writer, sheet_name='all_candidates')
if use_blast and blast_mirdeep_orig is not None:
    blast_mirdeep_orig.to_excel(writer, sheet_name='blast_miRdeep')
    blast_sRNAbench_orig.to_excel(writer, sheet_name='blast_sRNAbench')
writer.save()
