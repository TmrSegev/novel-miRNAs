from collections import defaultdict

from Bio import SeqIO
import re
import os
import configparser
import RNA
# import errors
import pandas as pd
from subprocess import Popen
from Bio import pairwise2
os.chdir(r"/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/")

dirpath = os.getcwd()
res = {}
def build_global_variables(): #modified
    # default parameters for the pipeline
    # mir35 - CACCGGG, let7 - GAGGUAG
    seed=short_window_size=long_window_size=max_energy=input_filter_parameters=organism_name_in_db=None
    path = ''
    windows = [65, 60, 55, 50, 45]
    param = 30
    # read settings file
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read(f'{path}settings.ini')

    if not settings.has_section('mode_1'):
        # print("mode 1 error: " + errors.errorMessage('section_not_exist'));
        print('OK continue') #modified
    # if settings.has_option('mode1_new', 'fasta_file_path'):
    #     fastaFile = settings.get('mode1_new', 'fasta_file_path')
    #     fastaFile.encode('ascii', 'ignore')
    # else:
    #     print("mode 1 error: " + errors.errorMessage('missing_input_file'));
    #
    if settings.has_option('mode_1', 'seed'):
        seed = settings.get('mode_1', 'seed')
        # seed = seed.encode('ascii', 'ignore')

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

# seed,short_window_size,long_window_size,max_energy,input_filter_parameters,organism_name_in_db=build_global_variables()
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
        ct_df.loc[index_i][4]) - 1  # file starts from index 1, so  decrease 1 for index to align with index 0
    repair_index_start_star = repair_index_end_mature

    index_i = start_mature
    repair_index_start_mature = 0
    decreased_start_mature = False
    while int(ct_df.loc[index_i][4]) == 0:
        index_i += 1
        repair_index_start_mature += 1
        decreased_start_mature = True
    # end_star = int(ct_df.loc[index_i][4])-1  # file starts from index 1, so  decrease 1 for index to align with index 0

    direct = False
    if int(ct_df.loc[start_mature - 2][4]) != 0:
        end_star_direct = int(ct_df.loc[start_mature - 2][4]) - 1
        direct = True
    else:
        end_star_undirect = int(
            ct_df.loc[index_i][4]) - 1  # file starts from index 1, so  decrease 1 for index to align with index 0

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
        while int(ct_df.loc[index_i+1][4]) == 0: #modified +1 because no index 0
            index_i += 1
            repair_index_start_mature += 1
            increased_start_mature = True
        end_hairpin = int(ct_df.loc[index_i+1][4]) - 1 #modified +1 because no index 0
        repair_index_end_star = repair_index_start_mature

        index_i = end_mature
        repair_index_end_mature = 0
        increased_end_mature = False
        while int(ct_df.loc[index_i][4]) == 0:
            index_i -= 1
            repair_index_end_mature += 1
            increased_end_mature = True
        # start_star = int(ct_df.loc[index_i][4])-1
    except Exception as e:
        print(e, index_i,repair_index_start_mature,repair_index_end_mature)
    direct = False
    if int(ct_df.loc[end_mature - 2][4]) != 0:
        start_star_direct = int(ct_df.loc[end_mature - 2][4]) - 1
        direct = True
    else:
        start_star_undirect = int(
            ct_df.loc[index_i][4]) - 1  # file starts from index 1, so  decrease 1 for index to align with index 0

    # fix problem that end of star at 5p is after the calculated end_hairpin
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
        if row[4] != '0': #modified
            if mature_df[0].iloc[0] <= int(row[4]) <= mature_df[0].iloc[len(mature_df) - 1]:  # check only connection outside mature #modified (int)
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
            # check only connection outside star
            if star_df[0].iloc[0] <= int(row[4]) <= star_df[0].iloc[
                len(star_df) - 1]:  # check only connection outside star
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

    # skip all first rows with zero connections - not a bulge
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
            mature_i = row[0] - 1  # start of bulge
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
                diff_bulge = abs(mature_bulge - star_bulge)
                if diff_bulge > max_bulge_symmetry:
                    max_bulge_symmetry = diff_bulge

                star_i = star_j

    return max_bulge_symmetry


def find_loop_size_3p(ct_df, start_mature):
    if int(ct_df.loc[start_mature - 1][4]) != 0:
        start_loop = int(ct_df.loc[start_mature - 1][4]) + 2  # Overhead
    else:
        # start_loop = int(ct_df.loc[start_mature + 1][4]) + 3
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
        end_loop = int(ct_df.loc[end_mature + 1][4]) + 2  # Overhead
    else:
        # start_loop = int(ct_df.loc[start_mature + 1][4]) + 3
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
    param,windows,settings,seed, short_window_size, long_window_size, max_energy, input_filter_parameters, organism_name_in_db = build_global_variables() #modified
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    # if settings.has_option('mode1_new', 'output_file_name'):
    #     output_file_name = settings.get('mode1_new', 'output_file_name')
    #     output_file_name = output_file_name.encode('ascii', 'ignore')
    #
    # if settings.has_option('mode1_new', 'output_path'):
    #     output_path = settings.get('mode1_new', 'output_path')
    #     output_path = output_path.encode('ascii', 'ignore')

    if settings.has_option('mode_1', 'fasta_file_path'):
        fastaFile = settings.get('mode_1', 'fasta_file_path')
        # fastaFile = fastaFile.encode('ascii', 'ignore')

    i = 1
    records_dict = SeqIO.index(fastaFile, "fasta")

    s = seed
    # s = s.replace("U", "T") #modified

    genome_size = 0

    for key, record in records_dict.items():
        fasta = record.seq._data.decode('utf-8')
        genome_size += len(fasta)
        reverse_complement_fasta = record.seq.reverse_complement_rna()._data.decode('utf-8')
        strand_plus = "+"
        strand_minus = "-"
        i = find_seed_in_fasta(s, fasta, i, strand_plus, key)
        i = find_seed_in_fasta(s, reverse_complement_fasta, i, strand_minus, key)
        # res.setdefault(key, {})['seq'] = fasta  #modified
    # print("candidates found in this phase (find_candidates_by_seed): ", str(len(res)))
    # print("genome_size: ", genome_size)

    # return res


# help function to create the candidates from fasta record
def find_seed_in_fasta(seed, fasta, i, strand, key):
    param,windows,settings,seed, short_window_size, long_window_size, max_energy, input_filter_parameters, organism_name_in_db = build_global_variables()
    len_fasta = len(fasta)
    for m in re.finditer(seed, fasta):
        # 5p
        start = m.start() - 1  # m.start() - short_window_size
        end = m.end() + short_window_size + long_window_size
        seq = fasta[start:end].replace('T', 'U')
        if len(seq) > 0:
            if strand == "-":
                temp = end
                end = len_fasta - start
                start = len_fasta - temp
            obj = {'Chr': str(key), 'Start_hairpin': str(start), 'End_hairpin': str(end), 'Strand': str(strand),
                   'Hairpin_seq': str(seq)}
            res[str(i)] = obj
            i = i + 1
        # 3p
        start = m.start() - long_window_size
        end = m.end() + short_window_size
        seq = fasta[start:end].replace('T', 'U')
        if len(seq) > 0:
            # res[str(i)] = {'chr': str(key), 'start': str(start), 'end': str(end), strand: str(strand)}
            if strand == "-":
                temp = end
                end = len_fasta - start
                start = len_fasta - temp
            obj = {'Chr': str(key), 'Start_hairpin': str(start), 'End_hairpin': str(end), 'Strand': str(strand),
                   'Hairpin_seq': str(seq)}
            res[str(i)] = obj
            i = i + 1
    return i


def fold_candidates():
    for key in res:
        fold = RNA.fold(res[key]['Seq'])
        res[key]['Fold'] = str(fold[0])


def fold_candidates_by_seq(seq):
    fold = RNA.fold(seq)
    return str(fold[0])


def find_seed(seed, seq):
    for m in re.finditer(seed, seq):
        start = m.start() + 1
        end = m.end()
        return start, end


def count_kmers(sequence, k, minmax):
    """Count kmer occurrences in a given sequence.

    Parameters
    ----------
    read : string
        A single DNA sequence.
    k : int
        The value of k for which to count kmers.

    Returns
    -------
    counts : dictionary, {'string': int}
        A dictionary of counts keyed by their individual kmers (strings
        of length k).

    """
    # Start with an empty dictionary
    # counts = {}
    counts = defaultdict(int)

    # Calculate how many kmers of length k there are
    num_kmers = len(sequence) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = sequence[i:i+k]
        ## Add the kmer to the dictionary if it's not there
        # if kmer not in counts:
        #     counts[kmer] = 0
        # Increment the count for this kmer
        counts[kmer] += 1
    # Return the final counts
    if minmax == "max":
        return max(counts.values())
    elif minmax == "min":
        return min(counts.values())
    else:
        raise Exception("Wrong minmax value")


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
        # if (index % 500) == 0:
        #     print(index)
        index += 1

        # if index == 50:
        #     break

        hairpin = value['Hairpin_seq']

        if true_mature:  # modified because seed appear twice
            mature = true_mature
        else:
            m = hairpin.find(seed)
            if m <= 0:
                del res[key]
                continue
            mature = hairpin[m - 1:m + 21]

        # if len(hairpin) <= len(mature): #modified
        #     del res[key]
        #     continue

        # check if 3p
        if check_if_mature_3p(hairpin, mature):
            mature_3p = True
            mature_5p = False
        else:
            mature_5p = True
            mature_3p = False

        finish_filter = False

        # for window in windows:
        #     if mature_3p:
        #         window_hairpin = hairpin[long_window_size - window:len(hairpin)]
        #     if mature_5p:
        #         window_hairpin = hairpin[:len(hairpin) - (long_window_size - window)]

            # if len(window_hairpin) <= len(mature):
            #     break

        # fold = fold_candidates_by_seq(window_hairpin)
        fold = fold_candidates_by_seq(hairpin)

        with open('ct.txt', 'w') as infile:
            infile.write('>' + key + '\n')
            infile.write(hairpin + '\n')
            infile.write(fold)

        # cmd = "RNAfold --noPS ct.txt | b2ct > ct_file.ct"
        cmd = "~/.conda/envs/my_env/bin/RNAfold --noPS ct.txt | ~/.conda/envs/my_env/bin/b2ct > ct_file.ct"  #TODO change env name
        p = Popen(cmd, shell=True)
        p.communicate()

        ct_df = pd.read_csv('ct_file.ct', delimiter='\s+', header=None, names=[0, 1, 2, 3, 4, 5])
        ct_df = ct_df.iloc[1:]
        ct_df.astype({4: 'int'}).dtypes

        if len(ct_df) == 0:
            continue

        # find indexes of the seed, mature, and hairpin
        if true_mature:
            start_mature = hairpin.find(true_mature) + 1
            start_seed, end_seed = start_mature + 1, start_mature + 8
            end_mature = min(start_mature + len(true_mature) - 1, len(ct_df))
            if true_mature == "UGAUAUGUUGUUUGAAUGCCCCU":
                print("UGAUAUGUUGUUUGAAUGCCCCU", start_mature, end_mature)
                print("start_mature:", start_mature, "start + len:", start_mature + len(true_mature), "len(ct_df:", len(ct_df))
        else:
            print("row 511, this should NOT be printed")
            #start_seed, end_seed = find_seed(seed, window_hairpin)
            start_seed, end_seed = find_seed(seed, hairpin)

            start_mature = start_seed - 1
            end_mature = min(end_seed + 14, len(ct_df))

        # minimum 16 base pairings in duplex
        mature_df = ct_df.loc[start_mature:end_mature]

        max_bulge_symmetry = find_max_bulge_symmetry(mature_df)

        # if max_bulge_symmetry > dict_filter_params['max_bulge_symmetry']:
        #     continue

        mature_numbers_of_connections, mature_max_bulge = mature_complimentarity(mature_df)
        mature_bp_ratio = mature_numbers_of_connections / float(len(mature))

        # if mature_bp_ratio < 0.68 or mature_max_bulge > dict_filter_params[  #modified
        #     'max_mature_bulge']:  # 15 connections, dict_filter_params['min_mature_bp_ratio']
        #     continue

        if mature_3p:
            # find indexes of start and end of the hairpin from ct file
            hairpin_boundries = ct_file_parser_3p(ct_df, start_mature, end_mature, param)
        if mature_5p:
            hairpin_boundries = ct_file_parser_5p(ct_df, start_mature, end_mature, param)

        # if hairpin_boundries['valid'] is False: #modified
        #     continue

        # cut the hairpin with the new indexes
       # cutted_hairpin = window_hairpin[hairpin_boundries['start_hairpin']:hairpin_boundries['end_hairpin'] + 1]
        #cutted_hairpin = hairpin[hairpin_boundries['start_hairpin']:hairpin_boundries['end_hairpin'] + 1]

        # if len(cutted_hairpin) < dict_filter_params['min_trimmed_hairpin_length']: #modified
        #     continue

        # if mature_3p:
        #     # find the indexes of the loop
        #     start_loop, end_loop = find_loop_size_3p(ct_df, start_mature)
        # if mature_5p:
        #     start_loop, end_loop = find_loop_size_5p(ct_df, end_mature)
        #
        # loop_size = end_loop - start_loop

        # if loop_size < dict_filter_params['min_loop_length'] or loop_size > dict_filter_params['max_loop_length']:
        #     # del res[key]
        #     continue

        # if mature_3p:
        #     start_star = hairpin_boundries['start_star']
        #     end_star = hairpin_boundries['end_star']
        #     star_length = end_star - start_star
        #
        # if mature_5p:
        #     start_star = hairpin_boundries['start_star']
        #     end_star = hairpin_boundries['end_star']
        #     star_length = end_star - start_star

        if true_star != "nan":
            start_star = hairpin.find(true_star)
            star_length = len(true_star)
            end_star = hairpin.find(true_star) + star_length
        else:
            start_star = hairpin_boundries['start_star']
            end_star = hairpin_boundries['end_star']
            star_length = end_star - start_star

        # if star_length < dict_filter_params['min_star_length'] or star_length > dict_filter_params[
        #     'max_star_length']:
        #     continue

        star_df = ct_df.loc[start_star + 1:end_star + 1]
        star_numbers_of_connections, star_max_bulge = star_complimentarity(star_df)
        star_bp_ratio = star_numbers_of_connections / float(len(star_df)) if len(star_df)>0 else 0

        # if star_bp_ratio < dict_filter_params['min_star_bp_ratio'] or star_max_bulge > dict_filter_params[
        #     'max_star_bulge']:  # 15 connections, dict_filter_params['min_mature_bp_ratio']
        #     continue
        if true_star != "nan":
            star = true_star
        else:
            star = hairpin[start_star:end_star + 1]
        # if mature_3p:
        #     star = window_hairpin[start_star:end_star+1]
        # if mature_5p:
        #     star = window_hairpin[start_star:end_star+1]

        ##########################################
        # fold = fold_candidates_by_seq(cutted_hairpin)
        #
        # with open('ct.txt', 'w') as infile:
        #     infile.write('>' + key + '\n')
        #     infile.write(cutted_hairpin + '\n')
        #     infile.write(fold)
        #
        # cmd = "RNAfold --noPS ct.txt | b2ct > ct_file.ct"
        # p = Popen(cmd, shell=True)
        # p.communicate()
        #
        # ct_df = pd.read_csv('ct_file.ct', delimiter='\s+', header=None, names=[0, 1, 2, 3, 4, 5])
        # ct_df = ct_df.iloc[1:]
        # ct_df.astype({4: 'int'}).dtypes
        #########################################

        mature = mature_df[1].str.cat()

        if mature_3p:
            # find the indexes of the loop
            start_loop = len(star)
            end_loop = hairpin.find(mature)
        if mature_5p:
            start_loop = len(mature)
            end_loop = hairpin.find(star)

        loop_size = end_loop - start_loop
        if hairpin == "guggccgcguggcucaauuggauagagcaccugacuacggaucaggagguugcagguucgaguccugcgguggucgau":
            print(start_loop)
            print(end_loop)
            print(loop_size)

        min_one_mer_mature = count_kmers(mature, 1, "min") / len(mature)
        min_one_mer_hairpin = count_kmers(hairpin, 1, "min") / len(hairpin)
        max_one_mer_mature = count_kmers(mature, 1, "max") / len(mature)
        max_two_mer_mature = count_kmers(mature, 2, "max") / len(mature)
        max_one_mer_hairpin = count_kmers(hairpin, 1, "max") / len(hairpin)
        max_two_mer_hairpin = count_kmers(hairpin, 2, "max") / len(hairpin)

        # loop_seq = get_loop(ct_df, start_loop, end_loop, loop_size)
        res[key]['Mature_connections'] = mature_numbers_of_connections
        res[key]['Mature_BP_ratio'] = '%.2f' % mature_bp_ratio
        res[key]['Mature_max_bulge'] = '%.2f' % mature_max_bulge
        res[key]['Loop_length'] = loop_size
        res[key]['Valid mir'] = hairpin_boundries['valid']
        res[key]['Fold'] = fold

        res[key]['Mature'] = mature
        res[key]['Mature_Length'] = len(mature) #modified

        if mature_3p:
            res[key]['3p/5p'] = '3p'
            # res[key]['Mature_5p'] = ''
        if mature_5p:
            res[key]['3p/5p'] = '5p'
            # res[key]['Mature_3p'] = ''
        res[key]['Hairpin_seq_trimmed'] = hairpin

        res[key]['Star'] = star
        res[key]['Start_star'] = start_star
        res[key]['End_star'] = end_star
        res[key]['Star_length'] = star_length
        res[key]['Star_connections'] = star_numbers_of_connections
        res[key]['Star_BP_ratio'] = '%.2f' % star_bp_ratio
        res[key]['Star_max_bulge'] = '%.2f' % star_max_bulge
        res[key]['Hairpin_seq_trimmed_length'] = len(hairpin)
        #res[key]['Window'] = window

        res[key]['Max_bulge_symmetry'] = max_bulge_symmetry
        res[key]['min_one_mer_mature'] = min_one_mer_mature
        res[key]['min_one_mer_hairpin'] = min_one_mer_hairpin
        res[key]['max_one_mer_mature'] = max_one_mer_mature
        res[key]['max_two_mer_mature'] = max_two_mer_mature
        res[key]['max_one_mer_hairpin'] = max_one_mer_hairpin
        res[key]['max_two_mer_hairpin'] = max_two_mer_hairpin

        # res[key]['start_mature'] = start_mature
        # res[key]['end_mature'] = end_mature

        # finish_filter = True
        # break

        # if finish_filter is False:
        #     del res[key]
        #     continue

    # print("candidates found in this phase (filter_candidates): " + str(len(res)))


# def reverse_fold(hairpin_fold):
#     hairpin_fold_r = hairpin_fold[::-1]
#     start_replace_from = hairpin_fold_r.count(')')
#     hairpin_fold_r = hairpin_fold_r.replace("(", ")")
#     hairpin_fold_r = hairpin_fold_r.replace(")", "(", start_replace_from)
#     return hairpin_fold_r


def check_if_mature_3p(seq, mirna):
    m = seq.find(mirna)
    if m == -1:
        return False
    before_mirna = m
    after_mirna = len(seq) - (m + len(mirna))
    # check if 3p or 5p
    if before_mirna > after_mirna:
        return True
    else:
        return False


def collect_db_data():
    # read settings file
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    if settings.has_option('mode_1', 'organism_name_in_db'):
        organism_name_in_db = settings.get('mode_1', 'organism_name_in_db')
        # organism_name_in_db = organism_name_in_db.encode('ascii', 'ignore')
    if settings.has_option('mode_1', 'db_file_name'):
        db_file_name = settings.get('mode_1', 'db_file_name')
        # db_file_name = db_file_name.encode('ascii', 'ignore')
    df = pd.read_csv(os.path.join(dirpath, db_file_name), delimiter=',')
    print("organism len in positive DB", len(df[df['Organism'].str.contains(organism_name_in_db)]))
    return df[df['Organism'].str.contains(organism_name_in_db)]


def list_db_hairpin_id(df):
    result = []
    for row_db in df.iterrows():
        hp_id_db = row_db[1]['Hairpin_name']  # check this is hairpin-id and not hairpin_id

        result.append(hp_id_db)
    return result


# def get_max_score(alignments):
#     max_score = 0
#     for a in alignments:
#         if max_score == 0:
#             max_score = a[2]
#         else:
#             if a[2] > max_score:
#                 max_score = a[2]
#     return max_score


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
                                                len(hairpin_res)) - i:  # we can have 3 mismatches - threshold
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
        found = found_in_DB(res[key]['Hairpin_seq_trimmed'], real_db_df, relevant_hp_id_db)
        # print('key : {}, found_in_DB : {}, found_in_DB_id : {}'.format(key,found['found'],found['hairpin_id']))
        # res[key]['found_in_DB'] = found['found']
        res[key]['Found_in_DB_id'] = found['Hairpin_name']


def write_final_results_file():
    # read settings file
    settings = configparser.ConfigParser()
    settings._interpolation = configparser.ExtendedInterpolation()
    settings.read('settings.ini')

    if settings.has_option('mode_1', 'output_file_name'):
        output_file_name = settings.get('mode_1', 'output_file_name')
        # output_file_name = output_file_name.encode('ascii', 'ignore')
        print("output_file_name in write_final_results_file : ", output_file_name)
    if settings.has_option('mode_1', 'output_path'):
        output_path = settings.get('mode_1', 'output_path')
        # output_path = output_path.encode('ascii', 'ignore')
        print("output_path in write_final_results_file : ", output_path)

    df = pd.DataFrame.from_dict(res, orient='index').fillna(0)
    df.to_csv(output_path + '/' + output_file_name + ".csv", index=True)
    print("write {}".format(output_path + '/' + output_file_name + ".csv"))


def add_fold_energy_to_candidates():
    for key in res:
        fold = RNA.fold(res[key]['Hairpin_seq_trimmed'])
        # res[key]['Final_fold_energy'] = str(fold[1])
        # res[key]['Fold_cut_tails'] = str(fold[0])
        res[key]['Energy'] = str(fold[1])
        res[key]['Fold'] = str(fold[0])


def filterResultsByMaxEnergy():
    if settings.has_option('mode_1', 'input_filter_parameters'):
        input_filter_parameters = settings.get('mode_1', 'input_filter_parameters')

    f = open(input_filter_parameters, 'r')
    data = f.read()
    f.close()
    dict_filter_params = eval(data)

    for key in list(res.keys()):
        if float(res[key]['Energy']) > dict_filter_params[
            'max_energy']:  # or float(res[key]['final_fold_energy']) < dict_filter_params['min_energy']
            del res[key]
    print("candidates found in this phase: {}".format(len(res)))


# def filterResultsByMinimumLinks(min_links):
#     i = 0
#     for key in res:
#         if res[key]['mirna_2'] != '':
#             if int(res[key]['links_num']) < min_links:
#                 res[key]['mirna_2'] = ''
#             else:
#                 i = i + 1
#     for key in res.keys():
#         if res[key]['mirna_2']=='':
#             del res[key]
#     print ("candidates found in this phase: " + str(i))

# html output structure varibles
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
        <h1><i class="fas fa-dna"></i> &nbsp; Candidates Results - """  """</h1>
        <br>
        <table class="table table-striped">
            <thead>
                <tr>
                    <th></th>
		    <th>index key</th>
                    <th>Mature</th>
                    <th>5p/3p</th>
                    <th>Star</th>
                    <th>Hairpin</th>
					<th>Chromosome</th>
                    <th>Strand</th>
                    <th>Start hairpin</th>
                    <th>End hairpin</th>
                    <th>Free energy</th>
                    <th>Loop size</th>
                    <th>Mature_BP_ratio</th>
                    <th>Star_BP_ratio</th>
                    <th>Found in DB</th>                   
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
        # read settings file
        settings = configparser.ConfigParser()
        settings._interpolation = configparser.ExtendedInterpolation()
        settings.read('settings.ini')

        if settings.has_option('mode_1', 'output_file_name'):
            output_file_name = settings.get('mode_1', 'output_file_name')
            # output_file_name = output_file_name.encode('ascii', 'ignore')
            print("output_file_name in create_html_file : ", output_file_name)

        if settings.has_option('mode_1', 'output_path'):
            output_path = settings.get('mode_1', 'output_path')
            # output_path = output_path.encode('ascii', 'ignore')
            print("output_path in create_html_file : ", output_path)

        write_file = open(output_path + '/' + output_file_name + ".html", "w")
        write_file.write(start_html)
        i = 1
        for key in res:
            write_file.write("<tr>")
            write_file.write("<td>" + str(i) + "</td>\n")
            # write_file.write("<td>" + str(i) + " (" + key + ")" + "</td>\n")  # we changed (added) this
            write_file.write("<td>" + key + "</td>\n")  # we changed (added) this
            write_file.write("<td>" + res[key]['Mature'] + "</td>\n")
            write_file.write("<td>" + res[key]['3p/5p'] + "</td>\n")
            write_file.write("<td>" + res[key]['Star'] + "</td>\n")
            # write_file.write("<td>" + res[key]['Mature_5p'] + "</td>\n")
            write_file.write(
                """<td> <a target="_blank" href="http://nibiru.tbi.univie.ac.at/forna/forna.html?id=fasta&file=%3Eheader\\n""" +
                res[key]['Hairpin_seq_trimmed'] + "\\n" + res[key][
                    'Fold'] + """ "> <i class="far fa-image"></i> </a> </td>\n""")
            write_file.write("<td>" + res[key]['Chr'] + "</td>\n")
            write_file.write("<td>" + res[key]['Strand'] + "</td>\n")
            write_file.write("<td>" + res[key]['Start_hairpin'] + "</td>\n")
            write_file.write("<td>" + res[key]['End_hairpin'] + "</td>\n")
            # write_file.write("<td>" + res[key]['mature_type'] + "</td>\n")
            write_file.write("<td>" + str(round(float(res[key]['Energy']), 2)) + "</td>\n")
            write_file.write("<td>" + str(res[key]['Loop_length']) + "</td>\n")
            write_file.write("<td>" + str(res[key]['Mature_BP_ratio']) + "</td>\n")
            write_file.write("<td>" + str(res[key]['Star_BP_ratio']) + "</td>\n")
            # svm = '';
            # if enable_SVM:
            #     svm = str(res[key]['svm_result'])
            # write_file.write("<td>" + svm + "</td>\n")
            write_file.write("<td>" + res[key]['Found_in_DB_id'] + "</td>\n")
            write_file.write("</tr>")
            i = i + 1
        write_file.write(end_html)
        write_file.close()
    finally:
        write_file.close()
        print("write {}".format(output_path + '/' + output_file_name + ".html"))


# start = time.clock()
def start_filtering(seq, true_mature=None, true_star=None):
    # settings = configparser.ConfigParser()
    # settings._interpolation = configparser.ExtendedInterpolation()
    # settings.read('settings.ini')
    # if settings.has_option('mode_1', 'fasta_file_path'):
    #     fastaFile = settings.get('mode_1', 'fasta_file_path')
    # records_dict = SeqIO.index(fastaFile, "fasta")
    # seq=seq
    # for key, record in records_dict.items():
    #     seq = record.seq._data.decode('utf-8')

    obj = {'Chr': 'Rom_gen', 'Start_hairpin':0, 'End_hairpin': len(seq), 'Strand': str(seq),
           'Hairpin_seq': str(seq)}
    res['new'] = obj
    # find_candidates_by_seed()
    filter_candidates(true_mature, true_star)
    return res
# print(start_filtering())
# print(len(res))
#
# add_fold_energy_to_candidates()
# print(len(res))
#
# print("filter results by maximum energy")
# filterResultsByMaxEnergy()
# print(len(res))
#
# add_found_in_DB_to_res()
# print(len(res))
#
# write_final_results_file()
# create_html_file()
#
# elapsed = (time.clock() - start)
# print("Program executed in " + str(elapsed))
# print("finish running mode 1")




## Count 5p and 3p from csv output files from mode 1
# # read settings file
# settings = configparser.ConfigParser()
# settings._interpolation = configparser.ExtendedInterpolation()
# settings.read('settings.ini')

# def counter_report(csv_file):
#     df = pd.read_csv(csv_file)
#     ser = df['3p/5p'].value_counts()
#     return {'3p': ser['3p'], '5p': ser['5p']}
#
#
# if settings.has_option('mode_4', 'dir_result_mode1_path'):
#     dir_result_mode1_path = settings.get('mode_4', 'dir_result_mode1_path')
#
# # go over all files in mode1 to get the 3p and 5p candidates
# for dirpath, dirs, files in os.walk(dir_result_mode1_path):
#     for fname in files:
#         file_name = os.path.join(dirpath,fname)
#         if file_name.endswith('.csv'):
#             # file_name = 'result_files_mode1_mir35/results_caenorhabditis_briggsae.csv'
#             print('file_name {} report counter {}'.format(file_name, counter_report(file_name)))
#             # print("*********************************")


#ACAAGAG {'Chr': 'ROM_gen', 'Start_hairpin': '34', 'End_hairpin': '121', 'Strand': '+', 'Hairpin_seq': 'AACAAGAGUGAGAUCAUUUUGAAAGCUGAUU', 'Mature_connections': 0, 'Mature_BP_ratio': '0.00', 'Mature_max_bulge': '0.00', 'Loop_length': -20, 'Fold': '....((((((....))))))...........', 'Mature': 'AACAAGAGUGAGAUCAUUUUGA', '3p/5p': '5p', 'Hairpin_seq_trimmed': 'AACAAGAGUGAGAUCAUUUUGAAAGC', 'Star': 'AGAGUGAGAUCAUUUUGAAAGC', 'Start_star': 4, 'End_star': 25, 'Star_length': 21, 'Star_connections': 0, 'Star_BP_ratio': '0.00', 'Star_max_bulge': '0.00', 'Hairpin_seq_trimmed_length': 26, 'Window': 65, 'Max_bulge_symmetry': 0}
#AACAAGA {'Chr': 'ROM_gen', 'Start_hairpin': '33', 'End_hairpin': '120', 'Strand': '+', 'Hairpin_seq': 'UAACAAGAGUGAGAUCAUUUUGAAAGCUGAUU', 'Mature_connections': 0, 'Mature_BP_ratio': '0.00', 'Mature_max_bulge': '0.00', 'Loop_length': -18, 'Fold': '.....((((((....))))))...........', 'Mature': 'UAACAAGAGUGAGAUCAUUUUG', '3p/5p': '5p', 'Hairpin_seq_trimmed': 'UAACAAGAGUGAGAUCAUUUUGAAAGCU', 'Star': 'GAGUGAGAUCAUUUUGAAAGCU', 'Start_star': 6, 'End_star': 27, 'Star_length': 21, 'Star_connections': 1, 'Star_BP_ratio': '0.05', 'Star_max_bulge': '4.00', 'Hairpin_seq_trimmed_length': 28, 'Window': 65, 'Max_bulge_symmetry': 0}
#