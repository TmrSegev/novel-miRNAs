import os
import random

import matplotlib.pyplot as plt
import pandas as pd
import Ziv_Git
from Bio import SeqIO
import sys
import numpy as np


def get_seq_data(path, start_end_mark=False):
    seq = {}
    for seq_record in SeqIO.parse(path, "fasta"):
        if species != "Hofstenia":
            # seq_id = seq_record.description.split(';')
            seq_id = seq_record.description
            # seq_id = seq_id[0]
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
    fasta_file_path = {precursors}
    output_file_name = results_hofstenia
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
                  'Hairpin_seq_trimmed_length': [], 'Max_bulge_symmetry': [], 'min_one_mer_mature': [], 'min_one_mer_hairpin': [], 'max_one_mer_mature': [], 'max_two_mer_mature': [], 'max_one_mer_hairpin': [], 'max_two_mer_hairpin': [], 'Valid mir': []}

def build_exception_dict():
    return {'Chr': -1, 'Start_hairpin': -1, 'End_hairpin': -1, 'Strand': -1, 'Hairpin_seq': -1,
                  'Mature_connections': -1, 'Mature_BP_ratio': -1, 'Mature_max_bulge': -1, 'Loop_length': -1,
                  'Fold': -1, 'Mature': -1,'Mature_Length': -1, '3p/5p': -1, 'Hairpin_seq_trimmed': -1, 'Star': -1, 'Start_star': -1,
                  'End_star': -1, 'Star_length': -1, 'Star_connections': -1, 'Star_BP_ratio': -1, 'Star_max_bulge': -1,
                  'Hairpin_seq_trimmed_length': -1, 'Max_bulge_symmetry': -1, 'min_one_mer_mature': -1, 'min_one_mer_hairpin': -1, 'max_one_mer_mature': -1, 'max_two_mer_mature': -1, 'max_one_mer_hairpin': -1, 'max_two_mer_hairpin': -1, 'Valid mir': False}
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

def plot_series(series, ticks):
    series.plot.hist()
    plt.title(series.name)
    print("name:", series.name, "min:", series.min(), "max:", series.max())
    plt.xticks(np.arange(series.min(), series.max() + ticks, ticks))
    plt.savefig("./figures/{}.png".format(series.name), dpi=300)
    plt.clf()

if __name__ == '__main__':
    precursors = None
    mature = None
    star = None
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
        elif arg == '--star':
            star = sys.argv[i + 1]
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
    star = get_seq_data(star, start_end_mark=False)
    if species == "miRGeneDB":
        all_remaining = pd.DataFrame()
    elif species != "Hofstenia":
        all_remaining = pd.read_excel(all_remaining_path, sheet_name="all_candidates")
    else:
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
            out_dict=Ziv_Git.start_filtering(seq, true_mature=mature[name], true_star=star[name])
            # if i == 678:
            #     print(out_dict)
            for k,v in out_dict['new'].items():
                mirdb_dict[k].append(v)
            #ziv_features = ziv_features.append(mirdb_df)
        except Exception as e:
            #print(e,name,seq,out_dict)
            print(e, name, seq)
            exception_dict = build_exception_dict()
            for k,v in exception_dict.items():
                mirdb_dict[k].append(v)
            # row = pd.DataFrame(columns=list(mirdb_dict.keys()))
            # row.loc[0] = ''
            continue
    mirdb_df = pd.DataFrame(mirdb_dict)
    all_remaining.reset_index(inplace=True, drop=True)
    mirdb_df.reset_index(inplace=True, drop=True)
    output = pd.concat([all_remaining, mirdb_df], axis=1)
    # output.to_csv("all_remaining_after_ziv_{}.csv".format(species), index=False)

    output = output.astype({'Mature_BP_ratio': 'float', 'Mature_max_bulge': 'float', 'Star_BP_ratio': 'float', 'Star_max_bulge': 'float'})
    writer = pd.ExcelWriter('all_remaining_after_ziv_{}.xlsx'.format(species))
    output.to_excel(writer, sheet_name='(A) Unfiltered)', index=False)
    sum_fc_thres_ok = output[output['sum_FC_m > thres'] == 1].copy()
    sum_fc_thres_ok.to_excel(writer, sheet_name='(B) sum_FC>100', index=False)
    no_novel451 = sum_fc_thres_ok[sum_fc_thres_ok['novel451'] == 0].copy()
    no_novel451.to_excel(writer, sheet_name='(C) Novel451', index=False)
    structural = no_novel451[no_novel451['Valid mir'] == True].copy()
    structural = structural[(structural["Mature_connections"] >= 14) & (structural["Mature_connections"] <= 22)]
    structural = structural[(structural["Mature_BP_ratio"] >= 0.6) & (structural["Mature_BP_ratio"] <= 0.96)]
    structural = structural[(structural["Mature_max_bulge"] <= 4)]
    structural = structural[(structural["Loop_length"] >= 10) & (structural["Loop_length"] <= 25)]
    structural = structural[(structural["Mature_Length"] >= 19) & (structural["Mature_Length"] <= 26)]
    structural = structural[(structural["Star_length"] >= 19) & (structural["Star_length"] <= 25)]
    structural = structural[(structural["Star_connections"] >= 14) & (structural["Star_connections"] <= 23)]
    structural = structural[(structural["Star_BP_ratio"] >= 0.6) & (structural["Star_BP_ratio"] <= 0.96)]
    structural = structural[structural["Star_max_bulge"] <= 4]
    structural = structural[structural["Hairpin_seq_trimmed_length"] >= 53]
    structural = structural[structural["Max_bulge_symmetry"] <= 3]
    structural = structural[structural["min_one_mer_hairpin"] >= 0.1]
    structural = structural[structural["max_one_mer_hairpin"] <= 0.45]
    structural.to_excel(writer, sheet_name='(D) Structural Features', index=False)
    writer.save()

    if species != "Hofstenia":
        if species == "Elegans":
            mirgenedb = output[(output['Description_mirgenedb'] != '.') & (output['Valid mir'] == True)]
        else:
            mirgenedb = output[output['Valid mir'] == True]
        plot_series(mirgenedb['Hairpin_seq_trimmed_length'], 5.0)
        plot_series(mirgenedb['Mature_connections'], 1.0)
        plot_series(mirgenedb['Mature_BP_ratio'].astype('float'), 0.05)
        plot_series(mirgenedb['Mature_max_bulge'].astype('float'), 1.0)
        plot_series(mirgenedb['Loop_length'], 2.0)
        plot_series(mirgenedb['Mature_Length'], 1.0)
        plot_series(mirgenedb['Star_length'], 1.0)
        plot_series(mirgenedb['Star_connections'], 1.0)
        plot_series(mirgenedb['Star_BP_ratio'].astype('float'), 0.05)
        plot_series(mirgenedb['Star_max_bulge'].astype('float'), 1.0)
        plot_series(mirgenedb['Max_bulge_symmetry'], 1.0)
        plot_series(mirgenedb['min_one_mer_hairpin'], 0.05)
        plot_series(mirgenedb['max_one_mer_hairpin'], 0.05)
        # mirgenedb['End_hairpin'].plot.hist()
        # plt.xticks(np.arange(mirgenedb['End_hairpin'].min(), mirgenedb['End_hairpin'].max() + 1, 5.0))
        # plt.savefig("./figures/hairpin_length.png", dpi=300)
        # plt.clf()
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


