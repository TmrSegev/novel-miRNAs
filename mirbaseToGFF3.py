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

def is_5p_3p(precursor, sub_seq):
    distance_to_start = precursor.index(sub_seq)
    distance_to_end = len(precursor) - distance_to_start - len(sub_seq)
    if distance_to_end < distance_to_start:
        return '3p'
    else:
        return '5p'


annotation = gffpd.read_gff3("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/cel_mirbase.gff3")
gff = annotation.df

hairpin = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/animalsHairpin.fa"),
                        'fasta')
hairpin_dict = SeqIO.to_dict(hairpin)

mature = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/animalsMature.fa"),
                        'fasta')
mature_dict = SeqIO.to_dict(mature)

gff['attributes'] = gff['attributes'].str.replace('|', ';', regex=False)
gff['name'] = gff['attributes'].str.split(';', expand=True)[2]
gff['name'] = gff['name'].str.split('=', expand=True)[1]

new_gff = pd.DataFrame(columns=gff.columns)
new_pre_gff = pd.DataFrame(columns=gff.columns)
flag_5p3p = ''
for index, row in gff.iterrows():
    if row['type'] == "miRNA_primary_transcript":
        hairpinSeq = str(hairpin_dict[row['name']].seq)

        # get two next rows
        next1_row = gff.loc[index + 1]
        next2_row = gff.loc[index + 2]

        # figure out which is 3p, which is 5p, and if there's one or two.
        if next1_row['type'] == "miRNA" and next2_row['type'] == "miRNA":
            next1_seq = str(mature_dict[next1_row['name']].seq)
            next2_seq = str(mature_dict[next2_row['name']].seq)
            next1_5p_3p = next1_row['name'][-2:]
            next2_5p_3p = next2_row['name'][-2:]
            if next1_5p_3p == '5p':
                seq5p = next1_seq
                seq3p = next2_seq
                start_5p = next1_row['start']
                end_5p = next1_row['end']
                start_3p = next2_row['start']
                end_3p = next2_row['end']
            elif next1_5p_3p == '3p':
                seq5p = next2_seq
                seq3p = next1_seq
                start_5p = next2_row['start']
                end_5p = next2_row['end']
                start_3p = next1_row['start']
                end_3p = next1_row['end']
            # is the strand + or -?
            # cut according to slides by isana and update start/end
            if row['strand'] == '+':
                hairpinSeq = hairpinSeq[(start_5p - row['start']):(end_3p - row['start'] + 1)]
                row['start'] = start_5p
                row['end'] = end_3p
            elif row['strand'] == '-':
                hairpinSeq = hairpinSeq[(row['end'] - end_5p):(row['end'] - start_3p) + 1]
                row['start'] = start_3p
                row['end'] = end_5p
        elif next2_row['type'] == "miRNA_primary_transcript":
            # no need to cut the precursor. Just assign 5p/3p
            next1_seq = str(mature_dict[next1_row['name']].seq)
            flag_5p3p = is_5p_3p(hairpinSeq, next1_seq)

        # append trimmed sequence to new gff
        row['attributes'] = row['attributes'] + ';' + hairpinSeq
        new_gff = new_gff.append(row)
        new_pre_gff = new_pre_gff.append(row)

    elif row['type'] == "miRNA":
        if flag_5p3p == '':
            flag_5p3p = row['name'][-2:]
        row['attributes'] = row['attributes'] + ';' + str(mature_dict[row['name']].seq) + ';' + flag_5p3p
        flag_5p3p = ''
        new_gff = new_gff.append(row)

new_gff = new_gff.drop(['name'], axis=1)
new_pre_gff = new_pre_gff.drop(['name'], axis=1)

annotation.df = new_gff
annotation.to_gff3('cel_mirbase_seq.gff3')

annotation.df = new_pre_gff
annotation.to_gff3('cel_mirbase_pre_only.gff3')