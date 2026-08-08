import pandas as pd
import gffpandas.gffpandas as gffpd
import sys

# -----GETTING INPUTS-----
intersections_table_path = None
gff_path = None
for i in range(1, len(sys.argv), 2):
    arg = sys.argv[i]
    if arg == '--intersections-table':
        intersections_table_path = sys.argv[i + 1]
    elif arg == '--gff':
        gff_path = sys.argv[i + 1]
    elif arg == '--help' or arg == '-h':
        print(f'Manual:\n'
              f' --intersections-table <path>: path to bedtools intersection .bed file, either -a mirdeep and -b mirdeep or -a sRNAbench and -b sRNAbench.\n'
              f' --gff <path>: path to pre only gff file of mirdeep/sRNAbench.\n'
              )
        sys.exit()

# bedtools -wao emits 9 (A) + 9 (B) + 1 (overlap bp). Extra cols are tolerated.
intersections = pd.read_csv(
    intersections_table_path,
    sep='\t',
    header=None,
    usecols=range(18),
    names=[
        'Chr_a', '.1', 'pre_miRNA_a', 'Start_a', 'End_a', '.2', 'Strand_a', '.3', 'Description_a',
        'Chr_b', '.4', 'pre_miRNA_b', 'Start_b', 'End_b', '.5', 'Strand_b', '.6', 'Description_b',
    ],
)
annotation = gffpd.read_gff3(gff_path)
gff = annotation.df

# Check if intersections table is empty
if intersections.empty:
    print("Warning: Intersections table is empty. No overlaps to process.")
    print("GFF file will remain unchanged.")
    sys.exit(0)

# Drop self-hits and non-overlaps (B filled with '.' / -1 under -wao/-loj).
intersections = intersections[
    (intersections['Description_a'] != intersections['Description_b'])
    & (intersections['Chr_b'].astype(str) != '.')
]

if intersections.empty:
    print(
        "No non-self overlaps at the requested -f threshold; "
        "GFF left unlabeled (this is OK when candidates do not overlap)."
    )
    sys.exit(0)

# Create rc_m columns and rc_s columns
intersections['rc_m_a'] = intersections['Description_a'].str.split(';', expand=True)[1]
intersections['rc_m_b'] = intersections['Description_b'].str.split(';', expand=True)[1]
intersections['rc_s_a'] = intersections['Description_a'].str.split(';', expand=True)[2]
intersections['rc_s_b'] = intersections['Description_b'].str.split(';', expand=True)[2]

overlaps = intersections[intersections['Strand_a'] == intersections['Strand_b']]
senseAnti = intersections[intersections['Strand_a'] != intersections['Strand_b']]
print(
    f"Labeling from {len(overlaps)} same-strand and "
    f"{len(senseAnti)} opposite-strand intersect rows."
)

# Marking as overlap
duplicatesCheck = []
for index, row in overlaps.iterrows():
    if row['Description_a'] not in duplicatesCheck:
        gff['attributes'] = gff['attributes'].where(row['Description_a'] != gff['attributes'], gff['attributes'] + ';overlap')
        duplicatesCheck.append(row['Description_a'])
    if row['Description_b'] not in duplicatesCheck:
        gff['attributes'] = gff['attributes'].where(row['Description_b'] != gff['attributes'], gff['attributes'] + ';overlap')
        duplicatesCheck.append(row['Description_b'])

# Marking as sense / antisense
duplicatesCheck = []
for index, row in senseAnti.iterrows():
    if row['Description_a'] not in duplicatesCheck:
        if row['rc_m_a'] > row['rc_m_b']:
            gff['attributes'] = gff['attributes'].where(row['Description_a'] != gff['attributes'], gff['attributes'] + ';sense')
            gff['attributes'] = gff['attributes'].where(row['Description_b'] != gff['attributes'], gff['attributes'] + ';antisense')
        elif row['rc_m_a'] < row['rc_m_b']:
            gff['attributes'] = gff['attributes'].where(row['Description_a'] != gff['attributes'], gff['attributes'] + ';antisense')
            gff['attributes'] = gff['attributes'].where(row['Description_b'] != gff['attributes'], gff['attributes'] + ';sense')
        elif row['rc_m_a'] == row['rc_m_b']:
            # If the mature counts are equal, compare the star counts:
            if row['rc_s_a'] >= row['rc_s_b']:
                gff['attributes'] = gff['attributes'].where(row['Description_a'] != gff['attributes'], gff['attributes'] + ';sense')
                gff['attributes'] = gff['attributes'].where(row['Description_b'] != gff['attributes'], gff['attributes'] + ';antisense')
            else:
                gff['attributes'] = gff['attributes'].where(row['Description_a'] != gff['attributes'], gff['attributes'] + ';antisense')
                gff['attributes'] = gff['attributes'].where(row['Description_b'] != gff['attributes'], gff['attributes'] + ';sense')

        duplicatesCheck.extend([row['Description_a'], row['Description_b']])


annotation.df = gff
annotation.to_gff3(gff_path)
