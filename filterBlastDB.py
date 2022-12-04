from Bio import SeqIO

nematodes = SeqIO.parse(open("/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/BLAST_DB/Caenorhabditis_AND_entry_typeSequence.fasta"),
                        'fasta')
for fasta in nematodes:
    name = fasta.id
    if ("miRNA" in name) or ("microRNA" in name):
        print(f">{name}".replace('-', '_'))
        print(str(fasta.seq))

