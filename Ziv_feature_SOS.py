import os
import random

import matplotlib.pyplot as plt
import pandas as pd
import Ziv_Git
from Bio import SeqIO
import sys
import numpy as np
import pprint


def get_seq_data(path, start_end_mark=False):
    seq = {}
    for seq_record in SeqIO.parse(path, "fasta"):
        seq_id = seq_record.description
        if start_end_mark:
            seq[seq_id] = ('S' + "".join(str(seq_record.seq)) + 'E')
        else:
            seq[seq_id] = ("".join(str(seq_record.seq)))
    return seq


def load_double_mature_names(double_mature_path):
    """
    Load candidate names from double_mature.fasta file.
    Returns a set of candidate names.
    """
    double_mature_names = set()
    try:
        for seq_record in SeqIO.parse(double_mature_path, "fasta"):
            # The candidate name is in the description/header
            double_mature_names.add(seq_record.description)
        print(f"Loaded {len(double_mature_names)} candidate names from double_mature.fasta")
    except FileNotFoundError:
        print(f"Warning: double_mature.fasta file not found at {double_mature_path}")
    except Exception as e:
        print(f"Error loading double_mature.fasta: {e}")
    return double_mature_names


def crete_fasta(name, seq):
    with open('fasta_example.fa', 'w') as f:
        f.write(f'> {name}\n'
                f'{seq}\n')


def create_setting_ini(seed):
    mode_1 = f"""[mode_1]
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
    keys = ['Chr', 'Start_hairpin', 'End_hairpin', 'Strand', 'Hairpin_seq',
            'Mature', 'Start_mature', 'End_mature', 'Mature_Length', 'Mature_connections', 'Mature_BP_ratio',
            'Mature_max_bulge', 'Loop_length', 'Fold', '3p/5p', 'Hairpin_seq_trimmed', 'Star', 'Start_star',
            'End_star', 'Star_length', 'Star_connections', 'Star_BP_ratio', 'Star_max_bulge',
            'Hairpin_seq_trimmed_length', 'Mature_max_bulge_asymmetry', 'Star_max_bulge_asymmetry',
            'min_one_mer_mature', 'min_one_mer_hairpin', 'max_one_mer_mature', 'max_two_mer_mature',
            'max_one_mer_hairpin', 'max_two_mer_hairpin', '5p_overhang', '3p_overhang', 'Valid mir'
    ]
    return {f'{key}_ziv': [] for key in keys}


def build_exception_dict():
    exception_values = {
        'Chr': -1, 'Start_hairpin': -1, 'End_hairpin': -1, 'Strand': -1, 'Hairpin_seq': -1,
        'Mature': -1, 'Start_mature': -1, 'End_mature': -1, 'Mature_Length': -1,
        'Mature_connections': -1, 'Mature_BP_ratio': -1, 'Mature_max_bulge': -1,
        'Loop_length': -1, 'Fold': -1, '3p/5p': -1, 'Hairpin_seq_trimmed': -1,
        'Star': -1, 'Start_star': -1, 'End_star': -1, 'Star_length': -1,
        'Star_connections': -1, 'Star_BP_ratio': -1, 'Star_max_bulge': -1,
        'Hairpin_seq_trimmed_length': -1, 'Mature_max_bulge_asymmetry': -1,
        'Star_max_bulge_asymmetry': -1, 'min_one_mer_mature': -1,
        'min_one_mer_hairpin': -1, 'max_one_mer_mature': -1, 'max_two_mer_mature': -1,
        'max_one_mer_hairpin': -1, 'max_two_mer_hairpin': -1,
        '5p_overhang': -1, '3p_overhang': -1, 'Valid mir': False
    }
    return {f'{k}_ziv': v for k, v in exception_values.items()}



def find_seed(name, seq):
    """
    Pull the 7-nt seed (positions 2–8) directly from the mature miRNA
    rather than searching inside the precursor.
    """
    mat = mature.get(name)
    if mat is None:
        raise KeyError(f"No mature sequence for '{name}'")
    if len(mat) < 8:
        raise ValueError(f"Mature sequence for '{name}' is only {len(mat)} nt long")
    return mat[1:8]


def find_gen_seed(seq):
    start_mature_seq = 'ZZZZZ'
    return seq[seq.index(start_mature_seq) + 5 + 1:seq.index(start_mature_seq) + 5 + 8]


def find_neg_seed(seq):
    r = random.randint(0, 1)
    if r:
        start_mature_seq = len(seq) - 22
    else:
        start_mature_seq = 0
    return seq[start_mature_seq + 1:start_mature_seq + 8]


def clean(seq):
    for char in ['DDDDD', 'FFFFF', 'ZZZZZ', 'BBBBB']:
        seq = seq.replace(char, "")
    return seq.split("\n")[0]


def plot_series(series, ticks, output_dir="./"):
    # Ensure figures directory exists
    figures_dir = output_dir + "figures/"
    os.makedirs(figures_dir, exist_ok=True)
    
    series.plot.hist()
    plt.title(series.name)
    print("name:", series.name, "min:", series.min(), "max:", series.max())
    plt.xticks(np.arange(series.min(), series.max() + ticks, ticks), rotation=45)
    mean = series.mean()
    std = series.std()
    plt.axvline(mean, color="red")
    plt.axvline(mean + std, color="yellow")
    plt.axvline(mean - std, color="yellow")
    plt.savefig(figures_dir + "{}_{}.png".format(species, series.name), dpi=300)
    plt.clf()


def manual_change(df, new_genome=False):
    """
    Manually corrects specific mismatched hairpin sequences in a pandas DataFrame.

    This function applies a correction logic to each row, identifying problematic
    entries by matching their full hairpin sequence and replacing it with the
    corrected version.

    All sequence comparisons and replacements are done using RNA sequences (with 'U' instead of 'T').

    Args:
        df (pd.DataFrame): A DataFrame where each row represents a miRNA and must
                           contain at least the 'hairpinSeq' column.
        new_genome (bool): If True, also applies corrections specific to the new genome.

    Returns:
        pd.DataFrame: The modified DataFrame with corrected hairpin sequences.
    """

    # --- Define the corrections here ---
    # The key is the ENTIRE incorrect hairpin sequence.
    # The value is the ENTIRE corrected hairpin sequence.
    # NOTE: We define these with 'T' for readability and convert to 'U' in the code.
    corrections_t = {
        # For new-mir-novel17_2
        "TAATGGATGATTACTATTATAAATAGGACAAAGAAGGTTAAATGGTACATAAACATATAATGTAAAGTTGTTTTATAAAGTTTAATGGATGATTACTATTATT":
            "TAATGGATGATTACTATTATAAATAGGACAAAGAAGGTTAAATGGTACATAAACATATAATGTAAAGTTGTTTTATAAAGTTTAATGGATGATTAATATTAGG",

        # For miR-335A
        "TAGTATGTATCATTGAAGAACTTTATAAGGATTTTCATTGATGTATGCTTGG":
            "TAGTATGTATCATTGAAGAACCTTTATAAGGATTTTCATTGATGTATGCTTGG",

        # For new-mir-novel20_2
        "TAGTGTGTGTCTTTGATGAACCTTGTTAAGGTTCTTAATTGATGTATGATTAGG":
            "TAGTGTGTGTCTTTGATGAACCTTGTTAAGGTTCTTAATTGATATATGATTAGG",

        # For miR-2839
        "TCAAACAGAAGTTTTATGCACCAGGTTTGAAGTCATGCTTAAGGCATCTTGTGCAGAATT":
            "TCAAACAGAAGTTTTATGCACCAGGTTTGAAGTCATGCTTAAGGGATCTTGTGCAGAATT"
    }

    # Corrections for new genome
    corrections_t_new_genome = {
        # For new-mir-novel80_1
        # Mismatch: ...ATAATC[T]GTAG... -> ...ATAATC[A]GTAG...
        "AGTTATTTAATAATCTGTAGAATAACACATTCTGCTATTGTTATGTAACTGC":
            "AGTTATTTAATAATCAGTAGAATAACACATTCTGCTATTGTTATGTAACTGC",

        # For new-mir-novel107_1
        # Mismatch: ...CATTGAT[G]TAT... -> ...CATTGAT[A]TAT...
        "TAGTGTTTATCTTTGAAGAATCGTAAAACGATTTTTCATTGATGTATGCTCGG":
            "TAGTGTTTATCTTTGAAGAATCGTAAAACGATTTTTCATTGATATATGCTCGG",

        # For new-mir-novel89
        # Mismatch: TAG[C]ATGTAT... -> TAG[G]ATGTAT...
        "TAGCATGTATCTTTGAAGAATCTTATATGGTTTTTCATTGATGTATGCTTGG":
            "TAGGATGTATCTTTGAAGAATCTTATATGGTTTTTCATTGATGTATGCTTGG",

        # For new-mir-novel112_2
        # Mismatch: ...CTACCCAC[A] -> ...CTACCCAC[T]
        "TGGGTAGTTCTCATGGCTGCCATCAGAAATCAAAGTAAACAAATTTTAAAATAGTCATGAGAACTACCCACA":
            "TGGGTAGTTCTCATGGCTGCCATCAGAAATCAAAGTAAACAAATTTTAAAATAGTCATGAGAACTACCCACT",

        # For new-mir-novel54_7
        # Mismatch: ...ACTAA[CT]TT... -> ...ACTAA[CC]TT...
        "TAGTGTGTAACTTTGACTAACTTTATATGATTTTTCATTGTTGTATACTTGG":
            "TAGTGTGTAACTTTGACTAACCTTATATGATTTTTCATTGTTGTATACTTGG",

        # For new-mir-novel60_3_3
        # Mismatch: ...ACGTCA[T]TTGC... -> ...ACGTCA[A]TTGC...
        "TAGCATTTCAGTTTTAGAGAGAGATTGTGTTCTCACTCTTAACGTCATTTGCTTGC":
            "TAGCATTTCAGTTTTAGAGAGAGATTGTGTTCTCACTCTTAACGTCAATTGCTTGC"
    }

    # Convert all correction sequences from DNA (T) to RNA (U)
    # The keys are the incorrect hairpins, values are the corrected ones.
    corrections_u = {
        key.replace('T', 'U'): val.replace('T', 'U')
        for key, val in corrections_t.items()
    }
    
    # If new_genome is True, merge in the new genome corrections
    if new_genome:
        corrections_u_new_genome = {
            key.replace('T', 'U'): val.replace('T', 'U')
            for key, val in corrections_t_new_genome.items()
        }
        corrections_u.update(corrections_u_new_genome)

    print(corrections_u)

    # --- FIX ---
    # Define a function to apply to each row of the DataFrame.
    def correct_hairpin(row):
        # Get the hairpin sequence from the row, remove spaces for matching.
        original_hairpin_u = row['hairpinSeq'].replace(" ", "")

        # Check if this hairpin is in our correction dictionary
        if original_hairpin_u in corrections_u:
            # If it is, return the corrected version
            print(f"--- Correcting hairpin for miRNA containing '{row.get('miRNA_ID', 'N/A')}' ---")
            return corrections_u[original_hairpin_u]
        else:
            # Otherwise, return the original sequence unchanged
            return row['hairpinSeq']

    print(f"Applying corrections to 'hairpinSeq' column...")

    # Use .apply() to execute the 'correct_hairpin' function on each row.
    # The result is a new Series which we assign back to the 'hairpinSeq' column.
    df['hairpinSeq'] = df.apply(correct_hairpin, axis=1)

    print("Corrections complete.")
    return df

if __name__ == '__main__':
    precursors = None
    mature = None
    star = None
    species = None
    all_remaining_path = None
    new_genome = False
    args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--new-genome':
            new_genome = True
            i += 1
            continue
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
                  f' --all-remaining <path>: path to the all remaining filtered csv file.\n'
                  f' --new-genome: use new genome folder structure for output files.\n')
            sys.exit()
        i += 2

    # Determine output directory based on new_genome flag
    if new_genome:
        output_dir = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/"
    else:
        output_dir = "./"

    precursors = get_seq_data(precursors, start_end_mark=False)
    mature = get_seq_data(mature, start_end_mark=False)
    star = get_seq_data(star, start_end_mark=False)
    if species == "miRGeneDB":
        all_remaining = pd.DataFrame()
        # Load double_mature names when processing mirgenedb
        double_mature_path = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/mirgenedb_data_v3/double_mature.fasta"
        double_mature_names = load_double_mature_names(double_mature_path)
    else:
        all_remaining = pd.read_excel(all_remaining_path, sheet_name="all_candidates")
        double_mature_names = set()

    gen_dict = build_dict()
    neg_dict = build_dict()
    mirdb_dict = build_dict()
    # Track candidate names for mirgenedb processing
    candidate_names_list = []
    for name, seq in precursors.items():
        try:
            seed = find_seed(name, seq)
        except Exception as e:
            print("Skipping", name, "–", e)
            continue
        print("seed is:", seed)
        create_setting_ini(seed)
        try:
            out_dict = Ziv_Git.start_filtering(seq, true_mature=mature[name], true_star=star[name])
            for k, v in out_dict['new'].items():
                mirdb_dict[k].append(v)
            # Track candidate name for successful processing
            if species == "miRGeneDB":
                candidate_names_list.append(name)
        except Exception as e:
            print("FAILED in start_filtering or append:", e)
            exception_dict = build_exception_dict()
            for k, v in exception_dict.items():
                mirdb_dict[k].append(v)
            # Track candidate name even for exceptions
            if species == "miRGeneDB":
                candidate_names_list.append(name)
            continue

    mirdb_df = pd.DataFrame(mirdb_dict)
    all_remaining.reset_index(inplace=True, drop=True)
    mirdb_df.reset_index(inplace=True, drop=True)
    output = pd.concat([all_remaining, mirdb_df], axis=1)

    # Add two_matures flag column for mirgenedb
    if species == "miRGeneDB":
        if len(candidate_names_list) == len(output):
            # Create flag column: 1 if candidate name is in double_mature, 0 otherwise
            output['two_matures'] = [1 if name in double_mature_names else 0 for name in candidate_names_list]
            print(f"Added two_matures flag: {sum(output['two_matures'])} candidates flagged as two matures")
            
            # Add fasta name and sequences from input fasta files
            output['fasta_name'] = candidate_names_list
            output['hairpin_seq_fasta'] = [precursors.get(name, '') for name in candidate_names_list]
            output['mature_seq_fasta'] = [mature.get(name, '') for name in candidate_names_list]
            output['star_seq_fasta'] = [star.get(name, '') for name in candidate_names_list]
            print(f"Added fasta name and sequence columns from input fasta files")
        else:
            print(f"Warning: Number of candidate names ({len(candidate_names_list)}) doesn't match output DataFrame length ({len(output)})")
            output['two_matures'] = 0
            # Still add columns even if lengths don't match (with empty/default values)
            # Pad or truncate candidate_names_list to match output length
            padded_names = candidate_names_list[:len(output)] if len(candidate_names_list) >= len(output) else candidate_names_list + [''] * (len(output) - len(candidate_names_list))
            output['fasta_name'] = padded_names
            output['hairpin_seq_fasta'] = [precursors.get(name, '') if name else '' for name in padded_names]
            output['mature_seq_fasta'] = [mature.get(name, '') if name else '' for name in padded_names]
            output['star_seq_fasta'] = [star.get(name, '') if name else '' for name in padded_names]
            print(f"Added fasta name and sequence columns from input fasta files (with padding/truncation due to length mismatch)")

    if species == "miRGeneDB" or species == "Hofstenia" or species == "Hofstenia_newGenome":
        output = output.astype(
            {'Mature_BP_ratio_ziv': 'float', 'Mature_max_bulge_ziv': 'float', 'Star_BP_ratio_ziv': 'float',
             'Star_max_bulge_ziv': 'float'})
        if species == "Hofstenia" or species == "Hofstenia_newGenome":
            output['Description'] = output[['Description_mirdeep', 'Description_sRNAbench']].astype(str).agg('__'.join, axis=1).str.replace(';', '|', regex=False).str.replace('ID=', '', regex=False).str.replace('.', '', regex=False)
            output = manual_change(output, new_genome=new_genome)
        writer = pd.ExcelWriter(output_dir + 'all_remaining_after_ziv_{}.xlsx'.format(species))
        output.to_excel(writer, sheet_name='(A) Unfiltered', index=False)
        writer.save()
    if species != "miRGeneDB" and species != "Hofstenia" and species != "Hofstenia_newGenome":
        output = output.astype(
            {'Mature_BP_ratio_ziv': 'float', 'Mature_max_bulge_ziv': 'float', 'Star_BP_ratio_ziv': 'float',
             'Star_max_bulge_ziv': 'float'})
        writer = pd.ExcelWriter(output_dir + 'all_remaining_after_ziv_{}.xlsx'.format(species))
        output.to_excel(writer, sheet_name='(A) Unfiltered', index=False)
        sum_fc_thres_ok = output[output['sum_FC_m > thres'] == 1].copy()
        sum_fc_thres_ok.to_excel(writer, sheet_name='(B) sum_FC>100', index=False)
        no_novel451 = sum_fc_thres_ok[sum_fc_thres_ok['novel451'] == 0].copy()
        no_novel451.to_excel(writer, sheet_name='(C) Novel451', index=False)
        structural = no_novel451[no_novel451['Valid mir_ziv'] == True].copy()
        structural = structural[
            (structural["Mature_connections_ziv"] >= 14) & (structural["Mature_connections_ziv"] <= 22)]
        structural = structural[
            (structural["Mature_BP_ratio_ziv"] >= 0.6) & (structural["Mature_BP_ratio_ziv"] <= 0.96)]
        structural = structural[(structural["Mature_max_bulge_ziv"] <= 4)]
        structural = structural[(structural["Loop_length_ziv"] >= 10) & (structural["Loop_length_ziv"] <= 25)]
        structural = structural[(structural["Mature_Length_ziv"] >= 19) & (structural["Mature_Length_ziv"] <= 26)]
        structural = structural[(structural["Star_length_ziv"] >= 19) & (structural["Star_length_ziv"] <= 25)]
        structural = structural[(structural["Star_connections_ziv"] >= 14) & (structural["Star_connections_ziv"] <= 23)]
        structural = structural[(structural["Star_BP_ratio_ziv"] >= 0.6) & (structural["Star_BP_ratio_ziv"] <= 0.96)]
        structural = structural[structural["Star_max_bulge_ziv"] <= 4]
        structural = structural[structural["Hairpin_seq_trimmed_length_ziv"] >= 53]
        structural = structural[structural["Max_bulge_symmetry_ziv"] <= 3]
        structural = structural[structural["min_one_mer_hairpin_ziv"] >= 0.1]
        structural = structural[structural["max_one_mer_hairpin_ziv"] <= 0.45]
        structural.to_excel(writer, sheet_name='(D) Structural Features', index=False)
        writer.save()

    if species != "miRGeneDB" and species != "Hofstenia" and species != "Hofstenia_newGenome":
        if species == "Elegans":
            mirgenedb = output[(output['Description_mirgenedb'] != '.') & (output['Valid mir_ziv'] == True)]
        else:
            mirgenedb = output[output['Valid mir_ziv'] == True]
        mirgenedb.to_csv("temp.csv", sep='\t', index=False)
        plot_series(mirgenedb['Hairpin_seq_trimmed_length_ziv'], 5.0, output_dir)
        plot_series(mirgenedb['Mature_connections_ziv'], 1.0, output_dir)
        plot_series(mirgenedb['Mature_BP_ratio_ziv'].astype('float'), 0.05, output_dir)
        plot_series(mirgenedb['Mature_max_bulge_ziv'].astype('float'), 1.0, output_dir)
        plot_series(mirgenedb['Loop_length_ziv'], 2.0, output_dir)
        plot_series(mirgenedb['Mature_Length_ziv'], 1.0, output_dir)
        plot_series(mirgenedb['Star_length_ziv'], 1.0, output_dir)
        plot_series(mirgenedb['Star_connections_ziv'], 1.0, output_dir)
        plot_series(mirgenedb['Star_BP_ratio_ziv'].astype('float'), 0.05, output_dir)
        plot_series(mirgenedb['Star_max_bulge_ziv'].astype('float'), 1.0, output_dir)
        plot_series(mirgenedb['Max_bulge_symmetry_ziv'], 1.0, output_dir)
        plot_series(mirgenedb['min_one_mer_hairpin_ziv'], 0.05, output_dir)
        plot_series(mirgenedb['max_one_mer_hairpin_ziv'], 0.05, output_dir)
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


