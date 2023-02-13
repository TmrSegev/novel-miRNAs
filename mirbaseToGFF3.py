# open gff file and fasta files
# iterate rows
# for every row in the gff, extract the name.
# if it's precursor check the precursor file, if mature check the mature fasta file
# look for the name in the fasta file.
# take the corresponding sequence
# save to gff using gffpandas

import pandas as pd
from Bio import SeqIO

gff = pd.read_csv("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/cel_mirbase.gff3", sep='\t', names=['Chr', '.1', 'pre/miRNA', 'Start', 'End', '.2', 'Strand', '.3', 'Description'],
                  skiprows=13)
print(gff)
hairpin = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/animalsHairpin.fa"),
                        'fasta')
hairpin_dict = SeqIO.to_dict(hairpin)
print(str(hairpin_dict['cel-let-7'].seq))
mature = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/mirbase_data/animalsMature.fa"),
                        'fasta')

gff['name'] = gff['Description'].str.split('|', expand=True)[2]
gff['name'] = gff['name'].str.split('=', expand=True)[1]

# for index, row in gff.iterrows():
    # if gff['pre/miRNA'] == "miRNA_primary_transcript":
    #     str
    # elif gff['pre/miRNA'] == "miRNA":

# for fasta in hairpin:
#     description = fasta.description
#     print(f">{description}")
#     print(str(fasta.seq))
