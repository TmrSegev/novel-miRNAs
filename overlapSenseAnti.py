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

intersections = pd.read_csv(intersections_table_path, sep='\t', names=['Chr_a', '.1', 'pre_miRNA_a', 'Start_a', 'End_a', '.2', 'Strand_a', '.3', 'Description_a', 'Chr_b', '.4', 'pre_miRNA_b', 'Start_b', 'End_b', '.5', 'Strand_b', '.6', 'Description_b'])
annotation = gffpd.read_gff3(gff_path)
gff = annotation.df

# Create rc_m columns and rc_s columns
intersections['rc_m_a'] = intersections['Description_a'].str.split(';', expand=True)[1]
intersections['rc_m_b'] = intersections['Description_b'].str.split(';', expand=True)[1]
intersections['rc_s_a'] = intersections['Description_a'].str.split(';', expand=True)[2]
intersections['rc_s_b'] = intersections['Description_b'].str.split(';', expand=True)[2]


intersections = intersections[intersections['Description_a'] != intersections['Description_b']]
overlaps = intersections[intersections['Strand_a'] == intersections['Strand_b']]
senseAnti = intersections[intersections['Strand_a'] != intersections['Strand_b']]

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
