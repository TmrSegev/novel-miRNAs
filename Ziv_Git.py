from collections import defaultdict

from Bio import SeqIO
import re
import os
import configparser
import RNA
import pandas as pd
from subprocess import Popen
from Bio import pairwise2
import numpy as np
os.chdir(r"/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/")
_ct_dump_count = 0

dirpath = os.getcwd()
res = {}
def build_global_variables():
    seed=short_window_size=long_window_size=max_energy=input_filter_parameters=organism_name_in_db=None
    path = ''
    windows = [65, 60, 55, 50, 45]
    param = 30
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read(f'{path}settings.ini')

    if not settings.has_section('mode_1'):
        print('OK continue')

    if settings.has_option('mode_1', 'seed'):
        seed = settings.get('mode_1', 'seed')

    if settings.has_option('mode_1', 'short_window_size'):
        short_window_size = int(settings.get('mode_1', 'short_window_size'))

    if settings.has_option('mode_1', 'long_window_size'):
        long_window_size = int(settings.get('mode_1', 'long_window_size'))

    if settings.has_option('mode_1', 'max_energy'):
        max_energy = float(settings.get('mode_1', 'max_energy'))

    if settings.has_option('mode_1', 'input_filter_parameters'):
        input_filter_parameters = settings.get('mode_1', 'input_filter_parameters')


    if settings.has_option('mode_1', 'organism_name_in_db'):
        organism_name_in_db = settings.get('mode_1', 'organism_name_in_db')
    return param,windows,settings,seed,short_window_size, long_window_size, max_energy, input_filter_parameters,organism_name_in_db

def ct_file_parser_3p(ct_df, start_mature, end_mature, param):
    index_i = end_mature
    repair_index_end_mature = 0
    decreased_end_mature = False
    valid = True
    while int(ct_df.loc[index_i][4]) == 0:
        index_i -= 1
        repair_index_end_mature += 1
        decreased_end_mature = True
    start_hairpin = int(
        ct_df.loc[index_i][4]) - 1
    repair_index_start_star = repair_index_end_mature

    index_i = start_mature
    repair_index_start_mature = 0
    decreased_start_mature = False
    while int(ct_df.loc[index_i][4]) == 0:
        index_i += 1
        repair_index_start_mature += 1
        decreased_start_mature = True

    direct = False
    if int(ct_df.loc[start_mature - 2][4]) != 0:
        end_star_direct = int(ct_df.loc[start_mature - 2][4]) - 1
        direct = True
    else:
        end_star_undirect = int(
            ct_df.loc[index_i][4]) - 1

    if start_hairpin > end_mature or start_hairpin > param:
        valid = False

    if decreased_end_mature:
        if decreased_start_mature:
            if direct:
                return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature), 'end_hairpin': end_mature - 1,
                        'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                        'end_star': end_star_direct, 'valid': valid}
            return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature), 'end_hairpin': end_mature - 1,
                    'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                    'end_star': end_star_undirect + repair_index_start_mature + 2, 'valid': valid}
        if direct:
            return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature), 'end_hairpin': end_mature - 1,
                    'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                    'end_star': end_star_direct, 'valid': valid}
        return {'start_hairpin': max(0, start_hairpin - repair_index_end_mature), 'end_hairpin': end_mature - 1,
                'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                'end_star': end_star_undirect + 2, 'valid': valid}
    if decreased_start_mature:
        if direct:
            return {'start_hairpin': start_hairpin, 'end_hairpin': end_mature - 1,
                    'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                    'end_star': end_star_direct, 'valid': valid}
        return {'start_hairpin': start_hairpin, 'end_hairpin': end_mature - 1,
                'start_star': max(0, start_hairpin - repair_index_start_star + 2),
                'end_star': end_star_undirect + 2, 'valid': valid}
    if direct:
        return {'start_hairpin': start_hairpin, 'end_hairpin': end_mature - 1, 'start_star': start_hairpin + 2,
                'end_star': end_star_direct, 'valid': valid}
    return {'start_hairpin': start_hairpin, 'end_hairpin': end_mature - 1, 'start_star': start_hairpin + 2,
            'end_star': end_star_undirect + 2, 'valid': valid}


def ct_file_parser_5p(ct_df, start_mature, end_mature, param):
    try:
        index_i = start_mature
        repair_index_start_mature = 0
        increased_start_mature = False
        valid = True
        while int(ct_df.loc[index_i+1][4]) == 0:
            index_i += 1
            repair_index_start_mature += 1
            increased_start_mature = True
        end_hairpin = int(ct_df.loc[index_i+1][4]) - 1
        repair_index_end_star = repair_index_start_mature

        index_i = end_mature
        repair_index_end_mature = 0
        increased_end_mature = False
        while int(ct_df.loc[index_i][4]) == 0:
            index_i -= 1
            repair_index_end_mature += 1
            increased_end_mature = True
    except Exception as e:
        print(e, index_i,repair_index_start_mature,repair_index_end_mature)
    direct = False
    if int(ct_df.loc[end_mature - 2][4]) != 0:
        start_star_direct = int(ct_df.loc[end_mature - 2][4]) - 1
        direct = True
    else:
        start_star_undirect = int(
            ct_df.loc[index_i][4]) - 1

    if min(len(ct_df), end_hairpin + repair_index_start_mature) < min(len(ct_df),
                                                                      end_hairpin + repair_index_end_star + 2):
        repair_index_start_mature = repair_index_end_star + 2

    if end_hairpin < start_mature or end_hairpin < param:
        valid = False

    if increased_start_mature:
        if increased_end_mature:
            if direct:
                return {'start_hairpin': start_mature - 1,
                        'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature),
                        'start_star': start_star_direct,
                        'end_star': min(len(ct_df), end_hairpin + repair_index_end_star + 2), 'valid': valid}
            return {'start_hairpin': start_mature - 1,
                    'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature),
                    'start_star': start_star_undirect - repair_index_end_mature + 2,
                    'end_star': min(len(ct_df), end_hairpin + repair_index_end_star + 2), 'valid': valid}
        if direct:
            return {'start_hairpin': start_mature - 1,
                    'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature),
                    'start_star': start_star_direct,
                    'end_star': min(len(ct_df), end_hairpin + repair_index_end_star + 2), 'valid': valid}
        return {'start_hairpin': start_mature - 1,
                'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature),
                'start_star': start_star_undirect + 2,
                'end_star': min(len(ct_df), end_hairpin + repair_index_end_star + 2), 'valid': valid}
    if increased_end_mature:
        if direct:
            return {'start_hairpin': start_mature - 1,
                    'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature),
                    'start_star': start_star_direct, 'end_star': end_hairpin + 2, 'valid': valid}
        return {'start_hairpin': start_mature - 1,
                'end_hairpin': min(len(ct_df), end_hairpin + repair_index_start_mature),
                'start_star': start_star_undirect - repair_index_end_mature + 2, 'end_star': end_hairpin + 2,
                'valid': valid}
    if direct:
        return {'start_hairpin': start_mature - 1, 'end_hairpin': end_hairpin, 'start_star': start_star_direct,
                'end_star': end_hairpin + 2, 'valid': valid}

    return {'start_hairpin': start_mature - 1, 'end_hairpin': end_hairpin, 'start_star': start_star_undirect + 2,
            'end_star': end_hairpin + 2, 'valid': valid}


def mature_complimentarity(mature_df):
    bulge_flag = False
    max_bulge = 0
    count_bulge = 0
    mature_connections = 0
    for row in mature_df.values:
        if row[4] != 0:
            if mature_df[0].iloc[0] <= int(row[4]) <= mature_df[0].iloc[len(mature_df) - 1]:
                continue
            mature_connections += 1

            if bulge_flag:
                if max_bulge < count_bulge:
                    max_bulge = count_bulge
                count_bulge = 0
                bulge_flag = False
        else:
            count_bulge += 1
            bulge_flag = True

    return mature_connections, max_bulge


def star_complimentarity(star_df):
    bulge_flag = False
    max_bulge = 0
    count_bulge = 0
    star_connections = 0
    for row in star_df.values:
        if int(row[4]) != 0:
            if star_df[0].iloc[0] <= int(row[4]) <= star_df[0].iloc[len(star_df) - 1]:
                continue
            star_connections += 1

            if bulge_flag:
                if max_bulge < count_bulge:
                    max_bulge = count_bulge
                count_bulge = 0
                bulge_flag = False
        else:
            count_bulge += 1
            bulge_flag = True

    return star_connections, max_bulge


def find_max_bulge_symmetry(mature_df):

    mature_i = 0
    mature_j = 0
    star_i = 0
    star_j = 0
    max_bulge_symmetry = 0
    seen_bulge = False

    counter = 0
    for row in mature_df.values:
        if int(row[4]) == 0:
            counter += 1
            continue
        else:
            break

    for index, row in enumerate(mature_df.values):
        if index < counter:
            continue
        if int(row[4]) == 0 and not seen_bulge:
            seen_bulge = True
            mature_i = row[0] - 1
        else:
            if int(row[4]) == 0:
                continue

            elif not seen_bulge:
                star_i = int(row[4])
                continue

            else:
                seen_bulge = False

                mature_j = row[0]
                star_j = int(row[4])

                mature_bulge = abs(mature_j - mature_i) - 1
                star_bulge = abs(star_j - star_i) - 1
                diff_bulge = mature_bulge - star_bulge
                if diff_bulge > max_bulge_symmetry:
                    max_bulge_symmetry = diff_bulge

                star_i = star_j

    return max_bulge_symmetry


def find_loop_size_3p(ct_df, start_mature):
    if int(ct_df.loc[start_mature - 1][4]) != 0:
        start_loop = int(ct_df.loc[start_mature - 1][4]) + 2
    else:
        index = start_mature + 1
        repair_index = 2
        while int(ct_df.loc[index][4]) == 0:
            index += 1
            repair_index += 1
        start_loop = int(ct_df.loc[index][4]) + repair_index

    end_loop = start_mature

    return start_loop, end_loop


def find_loop_size_5p(ct_df, end_mature):
    if int(ct_df.loc[end_mature + 1][4]) != 0:
        end_loop = int(ct_df.loc[end_mature + 1][4]) + 2
    else:
        index = end_mature - 1
        repair_index = 2
        while int(ct_df.loc[index][4]) == 0:
            index -= 1
            repair_index += 1
        end_loop = int(ct_df.loc[index][4]) - repair_index

    start_loop = end_mature

    return start_loop, end_loop


def get_loop(ct_df, start_loop, loop_size):
    loop = ''
    for i in range(loop_size):
        loop += ct_df.loc[start_loop + i][1]
    return loop


def find_candidates_by_seed():
    param,windows,settings,seed, short_window_size, long_window_size, max_energy, input_filter_parameters, organism_name_in_db = build_global_variables()
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    if settings.has_option('mode_1', 'fasta_file_path'):
        fastaFile = settings.get('mode_1', 'fasta_file_path')

    i = 1
    records_dict = SeqIO.index(fastaFile, "fasta")

    s = seed

    genome_size = 0

    for key, record in records_dict.items():
        fasta = record.seq._data.decode('utf-8')
        genome_size += len(fasta)
        reverse_complement_fasta = record.seq.reverse_complement_rna()._data.decode('utf-8')
        strand_plus = "+"
        strand_minus = "-"
        i = find_seed_in_fasta(s, fasta, i, strand_plus, key)
        i = find_seed_in_fasta(s, reverse_complement_fasta, i, strand_minus, key)


def find_seed_in_fasta(seed, fasta, i, strand, key):
    param,windows,settings,seed, short_window_size, long_window_size, max_energy, input_filter_parameters, organism_name_in_db = build_global_variables()
    len_fasta = len(fasta)
    for m in re.finditer(seed, fasta):
        # 5p
        start = m.start() - 1
        end = m.end() + short_window_size + long_window_size
        seq = fasta[start:end].replace('T', 'U')
        if len(seq) > 0:
            if strand == "-":
                temp = end
                end = len_fasta - start
                start = len_fasta - temp
            obj = {'Chr_ziv': str(key), 'Start_hairpin_ziv': str(start), 'End_hairpin_ziv': str(end), 'Strand_ziv': str(strand),
                   'Hairpin_seq_ziv': str(seq)}
            res[str(i)] = obj
            i = i + 1
        # 3p
        start = m.start() - long_window_size
        end = m.end() + short_window_size
        seq = fasta[start:end].replace('T', 'U')
        if len(seq) > 0:
            if strand == "-":
                temp = end
                end = len_fasta - start
                start = len_fasta - temp
            obj = {'Chr_ziv': str(key), 'Start_hairpin_ziv': str(start), 'End_hairpin_ziv': str(end), 'Strand_ziv': str(strand),
                   'Hairpin_seq_ziv': str(seq)}
            res[str(i)] = obj
            i = i + 1
    return i


def fold_candidates_by_seq(seq):
    fold = RNA.fold(seq)
    return str(fold[0])


def find_seed(seed, seq):
    for m in re.finditer(seed, seq):
        start = m.start() + 1
        end = m.end()
        return start, end


def count_kmers(sequence, k, minmax):
    counts = defaultdict(int)
    num_kmers = len(sequence) - k + 1
    for i in range(num_kmers):
        kmer = sequence[i:i+k]
        counts[kmer] += 1
    if minmax == "max":
        return max(counts.values())
    elif minmax == "min":
        return min(counts.values())
    else:
        raise Exception("Wrong minmax value")


def overhangs(mature_5p, ct_df, start_mature, end_mature, start_star, end_star):
    sm = start_mature - 1
    em = end_mature   - 1
    ss = start_star   - 1
    es = end_star     - 1

    if mature_5p:
        index = sm
        while ct_df.iloc[index, 4] == 0:
            index += 1
        index_on_star = ct_df.iloc[index, 4]
        overhang_3p = len(ct_df) - index_on_star - index

        index = ss
        shift = 0
        while ct_df.iloc[index, 4] == 0:
            index += 1
            shift += 1
        index_on_mature = ct_df.iloc[index, 4]
        overhang_5p = (em + 1) - index_on_mature - shift
    else:
        index = ss
        while ct_df.iloc[index, 4] == 0:
            index += 1
        index_on_mature = ct_df.iloc[index, 4]
        overhang_3p = len(ct_df) - index_on_mature - index

        index = sm
        shift = 0
        while ct_df.iloc[index, 4] == 0:
            index += 1
            shift += 1
        index_on_star = ct_df.iloc[index, 4]
        overhang_5p = (es + 1) - index_on_star - shift

    return overhang_5p, overhang_3p

def filter_candidates(true_mature=None, true_star=None):
    param,windows,settings,seed, short_window_size, long_window_size, max_energy, input_filter_parameters, organism_name_in_db = build_global_variables()
    if settings.has_option('mode_1', 'input_filter_parameters'):
        input_filter_parameters = settings.get('mode_1', 'input_filter_parameters')

    f = open(input_filter_parameters, 'r')
    data = f.read()
    f.close()
    dict_filter_params = eval(data)

    index = 0

    for key, value in list(res.items()):
        index += 1
        hairpin = value['Hairpin_seq_ziv']

        if true_mature:
            mature = true_mature
        else:
            m = hairpin.find(seed)
            if m <= 0:
                del res[key]
                continue
            mature = hairpin[m - 1:m + 21]

        if check_if_mature_3p(hairpin, mature):
            mature_3p = True
            mature_5p = False
        else:
            mature_5p = True
            mature_3p = False

        finish_filter = False

        fold = fold_candidates_by_seq(hairpin)

        with open('ct.txt', 'w') as infile:
            infile.write('>' + key + '\n')
            infile.write(hairpin + '\n')
            infile.write(fold)

        cmd = "~/.conda/envs/my_env/bin/RNAfold --noPS ct.txt | ~/.conda/envs/my_env/bin/b2ct > ct_file.ct"
        p = Popen(cmd, shell=True)
        p.communicate()

        ct_df = pd.read_csv('ct_file.ct', delimiter='\s+', header=None, names=[0, 1, 2, 3, 4, 5])
        ct_df = ct_df.iloc[1:]
        ct_df = ct_df.astype({4: 'int'})

        global _ct_dump_count
        if _ct_dump_count < 10:
            ct_df.to_csv(f'ct_debug_{_ct_dump_count}.tsv', sep='\t', index=False)
            print(f"[debug] wrote ct_debug_{_ct_dump_count}.tsv")
            _ct_dump_count += 1

        if len(ct_df) == 0:
            continue

        if true_mature:
            start_mature = hairpin.find(true_mature) + 1
            start_seed, end_seed = start_mature + 1, start_mature + 8
            end_mature = min(start_mature + len(true_mature) - 1, len(ct_df))
        else:
            start_seed, end_seed = find_seed(seed, hairpin)
            start_mature = start_seed - 1
            end_mature = min(end_seed + 14, len(ct_df))

        mature_df = ct_df.loc[start_mature:end_mature]

        mature_max_bulge_symmetry = find_max_bulge_symmetry(mature_df)
        mature_numbers_of_connections, mature_max_bulge = mature_complimentarity(mature_df)
        mature_bp_ratio = mature_numbers_of_connections / float(len(mature))

        if mature_3p:
            hairpin_boundries = ct_file_parser_3p(ct_df, start_mature, end_mature, param)
        if mature_5p:
            hairpin_boundries = ct_file_parser_5p(ct_df, start_mature, end_mature, param)

        if true_star and true_star != "nan":
            pos = hairpin.find(true_star)
            if pos == -1:
                raise ValueError(f"Star sequence {true_star} not found in hairpin")
            start_star = pos + 1
            star_length = len(true_star)
            end_star = start_star + star_length - 1
        else:
            start_star = hairpin_boundries['start_star']
            end_star = hairpin_boundries['end_star']
            star_length = end_star - start_star

        star_df = ct_df.loc[start_star:end_star + 1]
        star_max_bulge_symmetry = find_max_bulge_symmetry(star_df)
        star_numbers_of_connections, star_max_bulge = star_complimentarity(star_df)
        star_bp_ratio = star_numbers_of_connections / float(len(star_df)) if len(star_df)>0 else 0

        if true_star != "nan":
            star = true_star
        else:
            star = hairpin[start_star:end_star + 1]

        mature = mature_df[1].str.cat()

        if mature_3p:
            start_loop = len(star)
            end_loop = hairpin.find(mature)
        if mature_5p:
            start_loop = len(mature)
            end_loop = hairpin.find(star)

        loop_size = end_loop - start_loop

        min_one_mer_mature = count_kmers(mature, 1, "min") / len(mature)
        min_one_mer_hairpin = count_kmers(hairpin, 1, "min") / len(hairpin)
        max_one_mer_mature = count_kmers(mature, 1, "max") / len(mature)
        max_two_mer_mature = count_kmers(mature, 2, "max") / len(mature)
        max_one_mer_hairpin = count_kmers(hairpin, 1, "max") / len(hairpin)
        max_two_mer_hairpin = count_kmers(hairpin, 2, "max") / len(hairpin)

        overhang_5p, overhang_3p = overhangs(mature_5p, ct_df, start_mature, end_mature, start_star, end_star)

        res[key]['Mature_ziv'] = mature
        res[key]['Start_mature_ziv'] = start_mature
        res[key]['End_mature_ziv'] = end_mature
        res[key]['Mature_Length_ziv'] = len(mature)
        res[key]['Mature_connections_ziv'] = mature_numbers_of_connections
        res[key]['Mature_BP_ratio_ziv'] = '%.2f' % mature_bp_ratio
        res[key]['Mature_max_bulge_ziv'] = '%.2f' % mature_max_bulge
        res[key]['Loop_length_ziv'] = loop_size
        res[key]['Fold_ziv'] = fold


        if mature_3p:
            res[key]['3p/5p_ziv'] = '3p'
        if mature_5p:
            res[key]['3p/5p_ziv'] = '5p'
        res[key]['Hairpin_seq_trimmed_ziv'] = hairpin

        res[key]['Star_ziv'] = star
        res[key]['Start_star_ziv'] = start_star
        res[key]['End_star_ziv'] = end_star
        res[key]['Star_length_ziv'] = star_length
        res[key]['Star_connections_ziv'] = star_numbers_of_connections
        res[key]['Star_BP_ratio_ziv'] = '%.2f' % star_bp_ratio
        res[key]['Star_max_bulge_ziv'] = '%.2f' % star_max_bulge
        res[key]['Hairpin_seq_trimmed_length_ziv'] = len(hairpin)

        res[key]['Mature_max_bulge_asymmetry_ziv'] = mature_max_bulge_symmetry
        res[key]['Star_max_bulge_asymmetry_ziv'] = star_max_bulge_symmetry
        res[key]['min_one_mer_mature_ziv'] = min_one_mer_mature
        res[key]['min_one_mer_hairpin_ziv'] = min_one_mer_hairpin
        res[key]['max_one_mer_mature_ziv'] = max_one_mer_mature
        res[key]['max_two_mer_mature_ziv'] = max_two_mer_mature
        res[key]['max_one_mer_hairpin_ziv'] = max_one_mer_hairpin
        res[key]['max_two_mer_hairpin_ziv'] = max_two_mer_hairpin

        res[key]['5p_overhang_ziv'] = overhang_5p
        res[key]['3p_overhang_ziv'] = overhang_3p
        res[key]['Valid mir_ziv'] = hairpin_boundries['valid']


def check_if_mature_3p(seq, mirna):
    m = seq.find(mirna)
    if m == -1:
        return False
    before_mirna = m
    after_mirna = len(seq) - (m + len(mirna))
    if before_mirna > after_mirna:
        return True
    else:
        return False


def collect_db_data():
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    if settings.has_option('mode_1', 'organism_name_in_db'):
        organism_name_in_db = settings.get('mode_1', 'organism_name_in_db')
    if settings.has_option('mode_1', 'db_file_name'):
        db_file_name = settings.get('mode_1', 'db_file_name')
    df = pd.read_csv(os.path.join(dirpath, db_file_name), delimiter=',')
    print("organism len in positive DB", len(df[df['Organism'].str.contains(organism_name_in_db)]))
    return df[df['Organism'].str.contains(organism_name_in_db)]


def list_db_hairpin_id(df):
    result = []
    for row_db in df.iterrows():
        hp_id_db = row_db[1]['Hairpin_name']
        result.append(hp_id_db)
    return result


def found_in_DB(hairpin_res, real_db_df, relevant_hp_id_db):
    flag = False
    hp_id = ''
    mismatch_threshold = 3
    best_mismatches = mismatch_threshold
    for row_db in real_db_df.iterrows():
        hairpin_db = row_db[1]['Hairpin_seq_trimmed']
        alignments_best_score = pairwise2.align.globalms(hairpin_db, hairpin_res, 1, 0, 0, 0, penalize_end_gaps=False,
                                                         score_only=True)
        if alignments_best_score < min(len(hairpin_db), len(hairpin_res)) - mismatch_threshold:
            continue

        if row_db[1]['Hairpin_name'] in relevant_hp_id_db:
            for i in range(0, mismatch_threshold):
                if alignments_best_score >= min(len(hairpin_db),
                                                len(hairpin_res)) - i:
                    flag = True
                    min_mismatches = i
                    if min_mismatches < best_mismatches:
                        best_mismatches = min_mismatches
                        hp_id = row_db[1]['Hairpin_name']
                    break
    if flag:
        if hp_id in relevant_hp_id_db:
            relevant_hp_id_db.remove(hp_id)
    return {'found': flag, 'Hairpin_name': hp_id}


def add_found_in_DB_to_res():
    real_db_df = collect_db_data()
    relevant_hp_id_db = list_db_hairpin_id(real_db_df)
    for key in res:
        found = found_in_DB(res[key]['Hairpin_seq_trimmed_ziv'], real_db_df, relevant_hp_id_db)
        res[key]['Found_in_DB_id_ziv'] = found['Hairpin_name']


def write_final_results_file():
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    if settings.has_option('mode_1', 'output_file_name'):
        output_file_name = settings.get('mode_1', 'output_file_name')
        print("output_file_name in write_final_results_file : ", output_file_name)
    if settings.has_option('mode_1', 'output_path'):
        output_path = settings.get('mode_1', 'output_path')
        print("output_path in write_final_results_file : ", output_path)

    df = pd.DataFrame.from_dict(res, orient='index').fillna(0)
    df.to_csv(output_path + '/' + output_file_name + ".csv", index=True)
    print("write {}".format(output_path + '/' + output_file_name + ".csv"))


def add_fold_energy_to_candidates():
    for key in res:
        fold = RNA.fold(res[key]['Hairpin_seq_trimmed_ziv'])
        res[key]['Energy_ziv'] = str(fold[1])
        res[key]['Fold_ziv'] = str(fold[0])


def filterResultsByMaxEnergy():
    if settings.has_option('mode_1', 'input_filter_parameters'):
        input_filter_parameters = settings.get('mode_1', 'input_filter_parameters')

    f = open(input_filter_parameters, 'r')
    data = f.read()
    f.close()
    dict_filter_params = eval(data)

    for key in list(res.keys()):
        if float(res[key]['Energy_ziv']) > dict_filter_params[
            'max_energy']:
            del res[key]
    print("candidates found in this phase: {}".format(len(res)))


start_html = """<!DOCTYPE html>
<html lang="en">
<head>
    <style>
    table {
        text-align: center;
    }
    </style>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="ie=edge">
    <title>Document</title>
    <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.0.10/css/all.css" integrity="sha384-+d0P83n9kaQMCwj8F4RJB66tzIwOKmrdb46+porD/OvrJ+37WqIM7UoBtwHO6Nlg" crossorigin="anonymous">
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/css/bootstrap.min.css" integrity="sha384-9gVQ4dYFwwWSjIDZnLEWnxCjeSWFphJiwGPXr1jddIhOegiu1FwO5qRGvFXOdJZ4" crossorigin="anonymous">    
</head>
<body>
    <div class="container" style="max-width:1500px;">
        <br>
        <h1><i class="fas fa-dna"></i> &nbsp; Candidates Results - </h1>
        <br>
        <table class="table table-striped">
            <thead>
                <tr>
                    <th></th>
                    <th>index key</th>
                    <th>Mature_ziv</th>
                    <th>5p/3p_ziv</th>
                    <th>Star_ziv</th>
                    <th>Hairpin_ziv</th>
                    <th>Chromosome_ziv</th>
                    <th>Strand_ziv</th>
                    <th>Start hairpin_ziv</th>
                    <th>End hairpin_ziv</th>
                    <th>Free energy_ziv</th>
                    <th>Loop size_ziv</th>
                    <th>Mature_BP_ratio_ziv</th>
                    <th>Star_BP_ratio_ziv</th>
                    <th>Found in DB_ziv</th>                   
                </tr>
            </thead>
            <tbody>"""

end_html = """</tbody>
        </table>
    </div>
</body>
</html>"""


def create_html_file():
    try:
        settings = configparser.ConfigParser()
        settings._interpolation = configparser.ExtendedInterpolation()
        settings.read('settings.ini')

        if settings.has_option('mode_1', 'output_file_name'):
            output_file_name = settings.get('mode_1', 'output_file_name')
            print("output_file_name in create_html_file : ", output_file_name)

        if settings.has_option('mode_1', 'output_path'):
            output_path = settings.get('mode_1', 'output_path')
            print("output_path in create_html_file : ", output_path)

        write_file = open(output_path + '/' + output_file_name + ".html", "w")
        write_file.write(start_html)
        i = 1
        for key in res:
            write_file.write("<tr>")
            write_file.write("<td>" + str(i) + "</td>\n")
            write_file.write("<td>" + key + "</td>\n")
            write_file.write("<td>" + res[key]['Mature_ziv'] + "</td>\n")
            write_file.write("<td>" + res[key]['3p/5p_ziv'] + "</td>\n")
            write_file.write("<td>" + res[key]['Star_ziv'] + "</td>\n")
            write_file.write(
                """<td> <a target="_blank" href="http://nibiru.tbi.univie.ac.at/forna/forna.html?id=fasta&file=%3Eheader\\n""" +
                res[key]['Hairpin_seq_trimmed_ziv'] + "\\n" + res[key][
                    'Fold_ziv'] + """ "> <i class="far fa-image"></i> </a> </td>\n""")
            write_file.write("<td>" + res[key]['Chr_ziv'] + "</td>\n")
            write_file.write("<td>" + res[key]['Strand_ziv'] + "</td>\n")
            write_file.write("<td>" + res[key]['Start_hairpin_ziv'] + "</td>\n")
            write_file.write("<td>" + res[key]['End_hairpin_ziv'] + "</td>\n")
            write_file.write("<td>" + str(round(float(res[key]['Energy_ziv']), 2)) + "</td>\n")
            write_file.write("<td>" + str(res[key]['Loop_length_ziv']) + "</td>\n")
            write_file.write("<td>" + str(res[key]['Mature_BP_ratio_ziv']) + "</td>\n")
            write_file.write("<td>" + str(res[key]['Star_BP_ratio_ziv']) + "</td>\n")
            write_file.write("<td>" + res[key]['Found_in_DB_id_ziv'] + "</td>\n")
            write_file.write("</tr>")
            i = i + 1
        write_file.write(end_html)
        write_file.close()
    finally:
        write_file.close()
        print("write {}".format(output_path + '/' + output_file_name + ".html"))


def start_filtering(seq, true_mature=None, true_star=None):
    obj = {'Chr_ziv': 'Rom_gen', 'Start_hairpin_ziv':0, 'End_hairpin_ziv': len(seq), 'Strand_ziv': str(seq),
           'Hairpin_seq_ziv': str(seq)}
    res['new'] = obj
    filter_candidates(true_mature, true_star)
    return res