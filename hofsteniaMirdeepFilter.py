#!/usr/bin/python

import sys
import pandas as pd
import numpy as np
import io

pd.options.mode.chained_assignment = None

def filterInputs(inputs_arr, score_threshold, true_positive_threshold, mc_threshold, exclude_counts):
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
        to_delete = input[(input['mature read count'] < mc_threshold) & (input['star read count'] < mc_threshold)]
        to_delete['Removal Reason'] = 'max(mature read count, star read count) < {}'.format(mc_threshold)
        deleted_input = deleted_input.append(to_delete)
        input = input[(input['mature read count'] >= mc_threshold) | (input['star read count'] >= mc_threshold)]

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
        input.to_csv(f'remaining_file_{file_count - 1}.csv', sep='\t', index=False)
        filtered_inputs_arr.append(input)

        if deleted_input is not None:
            deleted_input.to_csv(f'removed_file_{file_count-1}.csv', sep='\t', index=False)
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

def run(inputs, threshold_tp, threshold_s, threshold_mc, exclude_c):
    f"""
    This function creates gff3 file from the inputs, This function will add to the miRNA ID -m / -s if its mature/star,
    -number which the number its the frequency of this seq among its type (mature/star)
    :param exclude_c: Int - exclude rows from the filtering when total counts higher than exclude_c.
    :param threshold_s: Float - threshold for score.
    :param threshold_tp: Float - threshold for the true positive estimate.
    :param inputs: Array - with DataFrames of miRDeep results.csv predict tables.
    :return: None, at the end of the function gff3 file will created.
    """
    filtered_input, removed_input = filterInputs(inputs, threshold_s, threshold_tp, threshold_mc, exclude_c)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    input = None
    csv_save = False
    threshold_tp = None
    threshold_s = None
    threshold_mc = None
    exclude_c = None
    args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        elif arg == '--filter-tp':
            threshold_tp = sys.argv[i + 1]
        elif arg == '--filter-s':
            threshold_s = sys.argv[i + 1]
        elif arg == '--filter-mc':
            threshold_mc = sys.argv[i + 1]
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
                  f' --filter-mc <float> : threshold for max counts, default: None.\n'
                  f' --exclude-c <int> : term to ignore the score filter threshold if total counts are higher, default: '
                  f'None.\n '
                  f' --csv-save : will save the inner tables of miRDeep2 output results as csv.\n'
                  f' --create-fasta <path>: create fasta file from the gff3 table.\n')
            sys.exit()
        i += 2

    if not input:
        raise ('Input path is required (-i <path>)')

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

    if threshold_mc is not None:
        threshold_mc = float(threshold_mc)

    if exclude_c is not None:
        exclude_c = int(exclude_c)

    run(inputs, threshold_tp, threshold_s, threshold_mc, exclude_c)
