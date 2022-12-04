from Bio import SeqIO

nematodes = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/BLAST_DB/Caenorhabditis_pre_miRNA.fasta"),
                        'fasta')
for fasta in nematodes:
    description = fasta.description
    print(f">{description}".replace(' ', '_'))
    print(str(fasta.seq))

