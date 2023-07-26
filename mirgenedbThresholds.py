from Bio import SeqIO
import sys
import pandas as pd


def separate_all_fasta(path):
    """
    Seperate sequeunces to 5 dictionaries, according to seq type.
    Save the 3 relevat dictionaries in to fasta files
    :param path:
    :return:
    """
    hairpin = {}
    mature = {}
    star = {}
    for seq_record in SeqIO.parse(path, "fasta"):
        seq_id = seq_record.description.split('_')[:-1][0]
        print(seq_id)
        suffix = seq_record.id.split('_')[-1]
        print(suffix)
        if suffix == "pre":
            hairpin[seq_id] = ("".join(str(seq_record.seq)))
        if '*' in suffix:
            star[seq_id] = ("".join(str(seq_record.seq)))
        elif suffix == "5p" or "3p":
            mature[seq_id] = ("".join(str(seq_record.seq)))
    return hairpin, mature, star


def save_dict_to_fasta(dictionary, output_file):
    with open(output_file, 'w') as file:
        for key, value in dictionary.items():
            file.write(f'>{key}\n{value}\n')
    print(f'Dictionary saved to {output_file}.')


hairpin, mature, star = separate_all_fasta("ALL.fas")
save_dict_to_fasta(hairpin, "ALL_mirgenedb_hairpin.fasta")
save_dict_to_fasta(mature, "ALL_mirgenedb_mature.fasta")
save_dict_to_fasta(star, "ALL_mirgenedb_star.fasta")
