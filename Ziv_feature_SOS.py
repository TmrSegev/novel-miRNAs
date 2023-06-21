import os
import random
import pandas as pd
import Ziv_Git
from Bio import SeqIO
import sys


def get_seq_data(path, start_end_mark=False):
    seq = {}
    for seq_record in SeqIO.parse(path, "fasta"):
        if species == "Elegans":
            seq_id = seq_record.id.split(';')
            print(seq_id)
            seq_id = seq_id[0]
        else:
            seq_id = seq_record.id.split('|')[:-1][0]
        if start_end_mark:
            seq[seq_id] = ('S' + "".join(str(seq_record.seq)) + 'E')
        else:
            seq[seq_id] = ("".join(str(seq_record.seq)))
    return seq

def crete_fasta(name,seq):
    with open('fasta_example.fa', 'w') as f:
        f.write(f'> {name}\n'
                f'{seq}\n')


def create_setting_ini(seed):
    #seed was CACCGGG
    mode_1=f"""[mode_1]
    seed = {seed}
    short_window_size = 14
    long_window_size = 65
    fasta_file_path = /sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_mirdeep_pre_only.fasta
    output_file_name = results_hofstenia_mirdeep
    output_path = .
    db_file_name = output_mir35_mode0_13_1_2023.csv
    organism_name_in_db = Caenorhabditis elegans
    input_filter_parameters = filter_parameters_mir35.txt
    """
    with open("settings.ini", "w") as f:
        f.write(mode_1)

def build_dict():
    return {'Chr': [], 'Start_hairpin': [], 'End_hairpin': [], 'Strand': [], 'Hairpin_seq': [],
                  'Mature_connections': [], 'Mature_BP_ratio': [], 'Mature_max_bulge': [], 'Loop_length': [],
                  'Fold': [], 'Mature': [],'Mature_Length': [], '3p/5p': [], 'Hairpin_seq_trimmed': [], 'Star': [], 'Start_star': [],
                  'End_star': [], 'Star_length': [], 'Star_connections': [], 'Star_BP_ratio': [], 'Star_max_bulge': [],
                  'Hairpin_seq_trimmed_length': [], 'Window': [], 'Max_bulge_symmetry': [], 'Valid mir': []}

def find_seed(name,seq):
    start_mature_inx = seq.index(mature[name])
    return seq[start_mature_inx + 1:start_mature_inx + 8]

def find_gen_seed(seq):
    start_mature_seq = 'ZZZZZ'
    return seq[seq.index(start_mature_seq) + 5 + 1:seq.index(start_mature_seq) + 5 + 8]


def find_neg_seed(seq):
    r = random.randint(0,1)
    if r:
        start_mature_seq = len(seq)-22
    else:
        start_mature_seq = 0
    return seq[start_mature_seq + 1:start_mature_seq + 8]

def clean(seq):
    for char in ['DDDDD','FFFFF','ZZZZZ','BBBBB']:
        seq = seq.replace(char,"")
    return seq.split("\n")[0]



if __name__ == '__main__':
    precursors = None
    mature = None
    species = None
    all_remaining_path = None
    args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--precursors':
            precursors = sys.argv[i + 1]
        elif arg == '--mature':
            mature = sys.argv[i + 1]
        elif arg == '--species':
            species = sys.argv[i + 1]
        elif arg == '--all-remaining':
            all_remaining_path = sys.argv[i + 1]
        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' --precursors <path> : fasta file path of precursors sequences.\n'
                  f' --mature <path> : fasta file path of mature sequences, with the same names as the precursors.\n'
                  f' --species <name>: name of the species.\n'
                  f' --all-remaining <path>: path to the all remaining filtered csv file.\n')
            sys.exit()
        i += 2

    precursors = get_seq_data(precursors, start_end_mark=False)
    mature = get_seq_data(mature, start_end_mark=False)
    all_remaining = pd.read_csv(all_remaining_path, sep='\t')
    #ziv_features = pd.DataFrame()
    # gen_seq = open("star_mature_generated_output.txt",).readlines()
    gen_dict = build_dict()
    neg_dict = build_dict()
    mirdb_dict = build_dict()
    for i,(name,seq) in enumerate(precursors.items()):
        seed = find_seed(name,seq)
        create_setting_ini(seed)
        try:
            out_dict=Ziv_Git.start_filtering(seq,true_mature=mature[name])
            if i == 678:
                print(out_dict)
            for k,v in out_dict['new'].items():
                mirdb_dict[k].append(v)
            #ziv_features = ziv_features.append(mirdb_df)
        except Exception as e:
            print(e,name,seq,out_dict)
            exception_dict = build_dict()
            for k,v in exception_dict.items():
                mirdb_dict[k].append(v)
            # row = pd.DataFrame(columns=list(mirdb_dict.keys()))
            # row.loc[0] = ''
            continue
    print("mir",len(mirdb_dict['Chr']),"\n")
    print(mirdb_dict)
    # #all_remaining = all_remaining.rename({"consensus precursor sequence":"Strand"}, axis=1)
    # print(all_remaining.info())
    mirdb_df = pd.DataFrame(mirdb_dict)
    print(mirdb_df.info())
    # print(output.info())
    print(all_remaining.info())
    #print(ziv_features.info())
    all_remaining.reset_index(inplace=True, drop=True)
    mirdb_df.reset_index(inplace=True, drop=True)
    output = pd.concat([all_remaining, mirdb_df], axis=1)
    output.to_csv("all_remaining_after_ziv_{}.csv".format(species))
    # for i,seq in enumerate(gen_seq):
    #     seed = find_gen_seed(seq)
    #     seq = clean(seq)
    #     create_setting_ini(seed)
    #     try:
    #         out_dict=Ziv_Git.start_filtering(seq)
    #         for k,v in out_dict['new'].items():
    #             gen_dict[k].append(v)
    #     except Exception as e:
    #         print(e,seq,out_dict)
    # print("gen",len(gen_dict['Chr']),"\n")
    # pd.DataFrame(gen_dict).to_csv("gen_df.csv")
    # fasta_neg_sequences = SeqIO.parse(open("e2e_negative_examples.fa"), 'fasta')
    # for i,fasta in enumerate(fasta_neg_sequences):
    #     if i>20000:
    #         break
    #     name, seq = fasta.id, str(fasta.seq)
    #     seed = find_neg_seed(seq)
    #     create_setting_ini(seed)
    #     try:
    #         out_dict=Ziv_Git.start_filtering(seq)
    #         for k,v in out_dict['new'].items():
    #             neg_dict[k].append(v)
    #     except Exception as e:
    #         print(i,e,name,seq,out_dict)
    # print("neg",len(neg_dict['Chr']))
    # pd.DataFrame(neg_dict).to_csv("neg_df.csv")


