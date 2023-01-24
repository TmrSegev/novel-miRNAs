#!/usr/bin/python

import sys
import pandas as pd
import numpy as np
import io

pd.options.mode.chained_assignment = None


def writeRemovedFasta(removed_input, seed_path, seed_file, fasta_path):
    open(fasta_path, 'w').close()
    removed_fasta_file = ''

    for input in removed_input:
        for index, row in input.iterrows():
            seq_id = getSeqId(row)
            star_seq = row['consensus star sequence']
            mature_seq = row['consensus mature sequence']
            hairpin = row['consensus precursor sequence']

            star_position = hairpin.index(star_seq)
            mature_position = hairpin.index(mature_seq)

            seq5p_id = seq_id + '|5p'
            seq3p_id = seq_id + '|3p'

            if star_position > mature_position:
                seq5p = row['consensus mature sequence']  # *
                seq3p = row['consensus star sequence']  # *
                seq5p_id += '|m'
                seq3p_id += '|s'

            else:
                seq5p = row['consensus star sequence']  # *
                seq3p = row['consensus mature sequence']  # *
                seq5p_id += '|s'
                seq3p_id += '|m'

            if seed_path:
                seq5p_seed = seq5p[1:8].upper()
                seq3p_seed = seq3p[1:8].upper()
                try:
                    seq5p_id += '|' + seed_file[seed_file['seed'] == seq5p_seed]["miRBase_name"].iloc[0]
                except:
                    seq5p_id += '|' + seq5p_seed

                try:
                    seq3p_id += '|' + seed_file[seed_file['seed'] == seq3p_seed]["miRBase_name"].iloc[0]
                except:
                    seq3p_id += '|' + seq3p_seed

            if seq5p != '-':
                removed_fasta_file += f'>{seq5p_id}\n{seq5p}\n'

            if seq3p != '-':
                removed_fasta_file += f'>{seq3p_id}\n{seq3p}\n'

            if len(removed_fasta_file) > 100000:
                with open(fasta_path, 'a+') as f:
                    f.write(removed_fasta_file)
                removed_fasta_file = ''

    with open(fasta_path, 'a+') as f:
        f.write(removed_fasta_file)


def filterInputs(inputs_arr, score_threshold, true_positive_threshold, exclude_counts):
    """
    This Function filtering the inputs Dataframe by threshold
    :param inputs_arr: inputs array
    :param score_threshold: Float - threshold for score.
    :param true_positive_threshold: Float - threshold for the true positive estimate.
    :return: filtered dataframes array
    """
    if exclude_counts is None:
        exclude_counts = float('inf')
    if score_threshold is None and true_positive_threshold is None:
        return inputs_arr

    filtered_inputs_arr = []
    removed_from_inputs_arr = []

    file_count = 1
    for input in inputs_arr:
        total_reads = len(input.index)
        total_reads_left = total_reads
        sys.stdout.write(f'Summary of Input File Number {file_count}:\n'
                         f'\tTotal Reads Before Filtering: {total_reads}\n')

        deleted_input = input[input['rfam alert'] != '-']
        deleted_input['Removal Reason'] = 'rfam alert'
        input = input[input['rfam alert'] == '-']
        for index, row in input.iterrows():
            with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_rRNA.fasta') as f:
                if row['consensus mature sequence'].upper() in f.read():
                    row['Removal Reason'] = 'rRNA'
                    deleted_input = deleted_input.append(row)
                    input.drop(index=index, inplace=True)
                    f.close()
            with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_snoRNA.fasta') as f:
                if row['consensus mature sequence'].upper() in f.read():
                    row['Removal Reason'] = 'snoRNA'
                    deleted_input = deleted_input.append(row)
                    input.drop(index=index, inplace=True)
                    f.close()
            with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_snRNA.fasta') as f:
                if row['consensus mature sequence'].upper() in f.read():
                    row['Removal Reason'] = 'snRNA'
                    deleted_input = deleted_input.append(row)
                    input.drop(index=index, inplace=True)
                    f.close()
            with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_tRNA.fasta') as f:
                if row['consensus mature sequence'].upper() in f.read():
                    row['Removal Reason'] = 'tRNA'
                    deleted_input = deleted_input.append(row)
                    input.drop(index=index, inplace=True)
                    f.close()
        total_reads_left = len(input.index)
        filtered_percent = "%.2f" % (((total_reads - total_reads_left) / total_reads) * 100)
        sys.stdout.write(
            f'\tTotal Reads Filtered by rfam alert and non-coding RNA database: {total_reads - total_reads_left} ({filtered_percent}%)\n')

        # filter the mature sequences that are already present with higher score.
        for index, row in input.iterrows():
            above_df = input.loc[:index, :].head(-1)
            if row['consensus mature sequence'] in above_df['consensus mature sequence'].values:
                if row['miRDeep2 score'] < score_threshold:  # if smaller than score threshold, remove the smaller one
                    row['Removal Reason'] = 'Has duplicate mature with higher score'
                    deleted_input = deleted_input.append(row)
                    input.drop(index=index, inplace=True)
            elif row['consensus mature sequence'] in above_df['consensus star sequence'].values:
                if row['miRDeep2 score'] < score_threshold:  # if smaller than score threshold, remove the smaller one
                    row['Removal Reason'] = 'Has duplicate star with higher score'
                    deleted_input = deleted_input.append(row)
                    input.drop(index=index, inplace=True)
        total_reads_left_after_duplicate_filtering = len(input.index)
        filtered_percent = "%.2f" % (
                ((total_reads_left - total_reads_left_after_duplicate_filtering) / total_reads_left) * 100)
        sys.stdout.write(
            f'\tTotal Reads Filtered because they have a duplicate with score >=10:'
            f' {total_reads_left - total_reads_left_after_duplicate_filtering} ({filtered_percent}%)\n')
        total_reads_left = total_reads_left_after_duplicate_filtering

        if score_threshold is not None:
            input['total read count'] = pd.to_numeric(input['total read count'])
            input['miRDeep2 score'] = pd.to_numeric(input['miRDeep2 score'])
            to_delete = input[((input['total read count'] < exclude_counts) | (input['star read count'] == 0)) & (input['miRDeep2 score'] < score_threshold)]
            to_delete['Removal Reason'] = 'Score < Threshold'
            deleted_input = deleted_input.append(to_delete)
            input = input[((input['total read count'] >= exclude_counts) & (input['star read count'] > 0)) | (input['miRDeep2 score'] >= score_threshold)]
            total_reads_left_after_score_filtering = len(input.index)
            filtered_percent = "%.2f" % (
                    ((total_reads_left - total_reads_left_after_score_filtering) / total_reads_left) * 100)
            sys.stdout.write(
                f'\tTotal Reads Filtered by Score:'
                f' {total_reads_left - total_reads_left_after_score_filtering} ({filtered_percent}%)\n')
            total_reads_left = total_reads_left_after_score_filtering
            #total_reads_left = len(input.index)
            #filtered_percent = "%.2f" % (((total_reads - total_reads_left) / total_reads) * 100)
           # sys.stdout.write(
            #    f'\tTotal Reads Filtered by Score: {total_reads - total_reads_left} ({filtered_percent}%)\n')

        if true_positive_threshold is not None:
            try:
                # Novel MicroRNA
                input['true positive probability'] = pd.to_numeric(input.apply(
                    lambda row: row['estimated probability that the miRNA candidate is a true positive'].split(' ')[0],
                    axis=1))
            except:
                # Known MicroRNA
                input['true positive probability'] = pd.to_numeric(input.apply(
                    lambda row: row['estimated probability that the miRNA is a true positive'].split(' ')[0],
                    axis=1))

            if deleted_input is not None:
                deleted_input = deleted_input.append(
                    input[(input['true positive probability'] < true_positive_threshold)])
            else:
                deleted_input = input[~(input['true positive probability'] >= true_positive_threshold)]
            input = input[input['true positive probability'] >= true_positive_threshold]
            total_reads_left_after_true_positive_filtering = len(input.index)
            filtered_percent = "%.2f" % (
                        ((total_reads_left - total_reads_left_after_true_positive_filtering) / total_reads_left) * 100)
            sys.stdout.write(
                f'\tTotal Reads Filtered by True-Positive Probability:'
                f' {total_reads_left - total_reads_left_after_true_positive_filtering} ({filtered_percent}%)\n')

            total_reads_left = total_reads_left_after_true_positive_filtering

        # filter mature / star
        to_delete = input[(input['mature read count'] < 100) & (input['star read count'] < 100)]
        to_delete['Removal Reason'] = 'max(mature read count, star read count) < 100'
        deleted_input = deleted_input.append(to_delete)
        input = input[(input['mature read count'] >= 100) | (input['star read count'] >= 100)]

        total_reads_left_after_mature_star_filtering = len(input.index)
        filtered_percent = "%.2f" % (
                ((total_reads_left - total_reads_left_after_mature_star_filtering) / total_reads_left) * 100)
        sys.stdout.write(
            f'\tTotal Reads Filtered because of low mature or star read count:'
            f' {total_reads_left - total_reads_left_after_mature_star_filtering} ({filtered_percent}%)\n')
        total_reads_left = total_reads_left_after_mature_star_filtering

        filtered_percent = "%.2f" % ((total_reads_left / total_reads) * 100)
        sys.stdout.write(f'\tTotal Reads Left After all filters: {total_reads_left} ({filtered_percent}%)\n')
        file_count += 1
        input.to_csv(f'remaining_file_{file_count - 1}.csv', sep='\t')
        filtered_inputs_arr.append(input)

        if deleted_input is not None:
            deleted_input.to_csv(f'removed_file_{file_count-1}.csv', sep='\t')
            removed_from_inputs_arr.append(deleted_input)
    return filtered_inputs_arr, removed_from_inputs_arr


def readMirbaseResults(input_path):
    """
    This function read the miRDeep2 results file (results_<date>.csv)
    and extract the inner tables of the predictions to DataFrames.
    :param input_path: String - the path of the results_<date>.csv.
    :return: Array - of Dataframes - which contain the prediction tables.
    """
    novel_string = ''
    mirbase_string = ''

    read_novel = False
    read_mirbase = False

    with open(input_path) as mirdeep_results:
        lines = mirdeep_results.readlines()

        for line in lines:
            if not read_mirbase:
                if read_novel and len(line) < 100:
                    read_novel = False
                if read_novel:
                    novel_string += line
                if 'novel miRNAs predicted by miRDeep2' in line:
                    read_novel = True

            if not read_novel:
                if read_mirbase and len(line) < 100:
                    read_mirbase = False
                if read_mirbase:
                    mirbase_string += line
                if 'mature miRBase miRNAs detected by miRDeep2' in line:
                    read_mirbase = True

    inputs = []
    try:
        novel_data = io.StringIO(novel_string)
        novel_df = pd.read_csv(novel_data, sep="\t")
        inputs.append(novel_df)
    except:
        pass

    try:
        mirbase_data = io.StringIO(mirbase_string)
        mirbase_df = pd.read_csv(mirbase_data, sep="\t")
        inputs.append(mirbase_df)
    except:
        pass

    return inputs


def getSeqId(row):
    try:
        seq_id = row['provisional id']  # *

    except:
        # That's mean that belong to known microRNA
        seq_id = row["mature miRBase miRNA"]
        seq_id = seq_id.replace('-3p', '')
        seq_id = seq_id.replace('-5p', '')

    return seq_id


def run(inputs, output, threshold_tp, threshold_s, exclude_c, fasta_path, seed_path):
    f"""
    This function creates gff3 file from the inputs, This function will add to the miRNA ID -m / -s if its mature/star,
    -number which the number its the frequency of this seq among its type (mature/star)
    :param seed_path: String - seed file path.
    :param fasta_path: String - fasta file will created in this path.
    :param exclude_c: Int - exclude rows from the filtering when total counts higher than exclude_c.
    :param threshold_s: Float - threshold for score.
    :param threshold_tp: Float - threshold for the true positive estimate.
    :param inputs: Array - with DataFrames of miRDeep results.csv predict tables.
    :param output: String - output path of the GFF formatted file.
    :return: None, at the end of the function gff3 file will created.
    """
    seed_file = None
    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)
    gff3_pre_only = pd.DataFrame(columns=gff3_columns)
    output_pre_only = "Elegans_mirdeep_pre_only.gff3"

    filtered_input, removed_input = filterInputs(inputs, threshold_s, threshold_tp, exclude_c)

    if seed_path is not None:
        seed_file = pd.read_csv(seed_path, sep="\t")

    if fasta_path is not None:
        fasta_file = ''
        open(fasta_path, 'w').close()

    intersection_index = -1  # Used later to intersect the table with miRdeep, blast and featurecounts results.
    for input in filtered_input:
        for index, row in input.iterrows():
            intersection_index += 1
            details = row['precursor coordinate']
            name = details.split(':')[0]
            positions = details.split(':')[1]
            strand = details.split(':')[2]
            seq_id = getSeqId(row)
            star_seq = row['consensus star sequence']
            mature_seq = row['consensus mature sequence']
            hairpin = row['consensus precursor sequence']
            rc_mature = row['mature read count']
            rc_star = row['star read count']

            star_position = hairpin.index(star_seq)
            mature_position = hairpin.index(mature_seq)

            seq5p_id = seq_id + '|5p'
            seq3p_id = seq_id + '|3p'

            if star_position > mature_position:
                seq5p = row['consensus mature sequence']  # *
                seq3p = row['consensus star sequence']  # *
                seq5p_freq = len(input[input['consensus mature sequence'] == seq5p])
                seq3p_freq = len(input[input['consensus star sequence'] == seq3p])
                seq5p_id += f'|m|{seq5p_freq}'
                seq3p_id += f'|s|{seq3p_freq}'
                mature_seq = 5

            else:
                seq5p = row['consensus star sequence']  # *
                seq3p = row['consensus mature sequence']  # *
                seq5p_freq = len(input[input['consensus star sequence'] == seq5p])
                seq3p_freq = len(input[input['consensus mature sequence'] == seq3p])
                seq5p_id += f'|s|{seq5p_freq}'
                seq3p_id += f'|m|{seq3p_freq}'
                mature_seq = 3

            seq5p_id += f'|index={intersection_index}'
            seq3p_id += f'|index={intersection_index}'

            if seed_path:
                if seq5p != '-':
                    seq5p_seed = seq5p[1:8].upper()
                    try:
                        seq5p_id += '|' + seed_file[seed_file['seed'] == seq5p_seed]["miRBase_name"].iloc[0]
                    except:
                        seq5p_id += '|' + seq5p_seed

                if seq3p != '-':
                    seq3p_seed = seq3p[1:8].upper()
                    try:
                        seq3p_id += '|' + seed_file[seed_file['seed'] == seq3p_seed]["miRBase_name"].iloc[0]
                    except:
                        seq3p_id += '|' + seq3p_seed

            if fasta_path is not None:
                if (seq5p != '-') & (mature_seq == 5):
                    fasta_file += f'>{seq5p_id}\n{seq5p}\n'

                if (seq3p != '-') & (mature_seq == 3):
                    fasta_file += f'>{seq3p_id}\n{seq3p}\n'

                if len(fasta_file) > 100000:
                    with open(fasta_path, 'a+') as f:
                        f.write(fasta_file)
                    fasta_file = ''

            start = int(positions.split('..')[0]) + 1
            end = int(positions.split('..')[1])
            if mature_seq == 5:
                seed = seq5p_id.split('|')[5]
                gff_row = [[f'{name}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={seq_id};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed}']]
            if mature_seq == 3:
                seed = seq3p_id.split('|')[5]
                gff_row = [[f'{name}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={seq_id};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed}']]
            gff3_pre_only = gff3_pre_only.append(gff_row)

            if strand == '+':
                if seq5p != '-':
                    offset5p = len(hairpin.split(seq5p)[0])
                    start5p = start + offset5p
                    end5p = start + offset5p + len(seq5p) - 1
                    gff_row.append([name, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={seq5p_id}'])

                if seq3p != '-':
                    offset3p = len(hairpin.split(seq3p)[0])
                    start3p = start + offset3p
                    end3p = start + offset3p + len(seq3p) - 1
                    gff_row.append([name, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={seq3p_id}'])

            else:
                if seq5p != '-':
                    offset5p = len(hairpin.split(seq5p)[0])
                    end5p = end - offset5p
                    start5p = end - offset5p - len(seq5p) + 1
                    gff_row.append([name, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={seq5p_id}'])

                if seq3p != '-':
                    offset3p = len(hairpin.split(seq3p)[0])
                    end3p = end - offset3p
                    start3p = end - offset3p - len(seq3p) + 1
                    gff_row.append([name, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={seq3p_id}'])

            miRNAs = pd.DataFrame(gff_row, columns=gff3_columns)

            gff3 = gff3.append(miRNAs)

    with open(output, 'w') as file:
        file.write(version)

    with open(output_pre_only, 'w') as file:
        file.write(version)

    if fasta_path is not None:
        with open(fasta_path, 'a+') as f:
            f.write(fasta_file)

        writeRemovedFasta(removed_input, seed_path, seed_file, 'removed_by_filter_' + fasta_path)

    gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')
    gff3_pre_only.to_csv(output_pre_only, index=False, header=False, mode="a", sep='\t')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    input = None
    output = None
    csv_save = False
    threshold_tp = None
    threshold_s = None
    fasta_path = None
    exclude_c = None
    seed_path = None
    args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        elif arg == '-o':
            output = sys.argv[i + 1]
        elif arg == '-seed':
            seed_path = sys.argv[i + 1]
        elif arg == '--filter-tp':
            threshold_tp = sys.argv[i + 1]
        elif arg == '--filter-s':
            threshold_s = sys.argv[i + 1]
        elif arg == '--create-fasta':
            fasta_path = sys.argv[i + 1]
        elif arg == '--exclude-c':
            exclude_c = sys.argv[i + 1]
        elif arg == '--csv-save':
            csv_save = True
            i += 1
            continue

        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i <path> : miRDeep2 prediction output path, like result_08_10_2021_t_09_57_05\n'
                  f' -o <path> : output path.\n'
                  f' -seed <path> : classify the reads by seed file, should be separated by tab with columns'
                  f' [miRBase_name, seed], default: None.\n'
                  f' --filter-tp <float> : threshold for the true positive estimate, any value between 0 - 100, '
                  f'default: None.\n '
                  f' --filter-s <float> : threshold for score, default: None.\n'
                  f' --exclude-c <int> : term to ignore the score filter threshold if total counts are higher, default: '
                  f'None.\n '
                  f' --csv-save : will save the inner tables of miRDeep2 output results as csv.\n'
                  f' --create-fasta <path>: create fasta file from the gff3 table.\n')
            sys.exit()
        i += 2

    if not input:
        raise ('Input path is required (-i <path>)')
    if not output:
        raise ('Output path is required (-o <path>)')

    inputs = readMirbaseResults(input)
    if csv_save is not None:
        count = 1
        for input in inputs:
            input.to_csv(f'table{count}.csv', sep='\t')
            count += 1

    if threshold_tp is not None:
        threshold_tp = float(threshold_tp)

    if threshold_s is not None:
        threshold_s = float(threshold_s)

    if exclude_c is not None:
        exclude_c = int(exclude_c)

    run(inputs, output, threshold_tp, threshold_s, exclude_c, fasta_path, seed_path)
