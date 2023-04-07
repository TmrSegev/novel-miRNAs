# open gff file and fasta files
# iterate rows
# for every row in the gff, extract the name.
# if it's precursor check the precursor file, if mature check the mature fasta file
# look for the name in the fasta file.
# take the corresponding sequence
# save to gff using gffpandas

import pandas as pd
from Bio import SeqIO
import gffpandas.gffpandas as gffpd

annotation = gffpd.read_gff3("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/cel_mirbase.gff3")
gff = annotation.df

hairpin = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/animalsHairpin.fa"),
                        'fasta')
hairpin_dict = SeqIO.to_dict(hairpin)

mature = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/animalsMature.fa"),
                        'fasta')
mature_dict = SeqIO.to_dict(mature)

gff['name'] = gff['attributes'].str.split('|', expand=True)[2]
gff['name'] = gff['name'].str.split('=', expand=True)[1]

new_gff = pd.DataFrame(columns=gff.columns)
for index, row in gff.iterrows():
    if row['type'] == "miRNA_primary_transcript":
        row['attributes'] = row['attributes'] + '|' + str(hairpin_dict[row['name']].seq)
        new_gff = new_gff.append(row)
    elif row['type'] == "miRNA":
        row['attributes'] = row['attributes'] + '|' + str(mature_dict[row['name']].seq)
        new_gff = new_gff.append(row)

new_gff = new_gff.drop(['name'], axis=1)

annotation.df = new_gff
annotation.to_gff3('cel_mirbase_seq.gff3')
