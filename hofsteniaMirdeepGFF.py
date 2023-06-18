#!/usr/bin/python

import sys
import pandas as pd
import numpy as np
import io


def getSeqId(row):
    try:
        seq_id = row['provisional id']  # *

    except:
        # That's mean that belong to known microRNA
        seq_id = row["mature miRBase miRNA"]
        seq_id = seq_id.replace('-3p', '')
        seq_id = seq_id.replace('-5p', '')

    return seq_id


def run(output, fasta_path, seed_path):
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
    output_pre_only = "{}_mirdeep_pre_only.gff3".format(species)
    fasta_pre_only_path = fasta_path.split('.fasta')[0]
    fasta_pre_only_path += "_pre_only.fasta"

    filtered_input = []
    for i in range(1, 3):
        # Uniting all remaining files
        table = None
        folders = ["EC1", "EC2", "EC3", "GA1", "GA2", "GA3", "DI1", "DI2", "DI3", "PDi1", "PDi2", "PDi3", "PDii1", "PDii2", "PDii3", "PL1", "PL2", "PL3", "PH1", "PH2", "PH3", "HL1", "HL2", "HL3", "IST1", "IST2", "IST3", "AMP1", "AMP2", "AMP3", "SMA1", "SMA2", "SMA3"]
        for folder in folders:
            to_add = pd.read_csv("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/mirdeep_out/" + folder + "/remaining_file_" + str(i) + ".csv", sep='\t')
            if table is None:
                table = to_add
            else:
                table = pd.concat([table, to_add], ignore_index=True)

        # Filtering by coordinates
        table = table.sort_values(['precursor coordinate'])

        table['chr'] = table['precursor coordinate'].str.split(':', expand=True)[0]
        table['positions'] = table['precursor coordinate'].str.split(':', expand=True)[1]
        table['positions'] = table['positions'].astype(str)
        table['start'] = table['positions'].str.split('.', expand=True)[0]
        table['start'] = table['start'].astype(int)
        table['end'] = table['positions'].str.split('.', expand=True)[2]
        table['end'] = table['end'].astype(int)
        table['overlaps'] = np.zeros(len(table))
        no_overlaps = pd.DataFrame(columns=table.columns)

        for index, row in table.iterrows():
            if index in table.index:
                table['distance'] = (row['end'] - table['start']) / (row['end'] - row['start'])
                overlaps = table[(table['distance'] >= 0.6) & (table['distance'] <= 1)].tail(-1)
                overlaps = overlaps[overlaps['chr'] == row['chr']]
                table.loc[index, 'overlaps'] = len(overlaps)
                if len(overlaps) == 0:
                    no_overlaps = no_overlaps.append(row)
                    table = table.drop(index)
                else:
                    table = table.drop(overlaps.index)
        print(table['overlaps'].value_counts().sort_index(ascending=False))
        filtered_input.append(table)
        no_overlaps.to_csv('removed_mirdeep_{}_no_overlaps.csv'.format(i), sep='\t')
        table = table.rename({"tag id":"provisional id", "estimated probability that the miRNA is a true positive":"estimated probability that the miRNA candidate is a true positive"}, axis=1)
    all_remaining = pd.concat([filtered_input[0], table], ignore_index=True)
    all_remaining.to_csv('mirdeep_all_remaining_filtered.csv', sep='\t', index=False)

    if seed_path is not None:
        seed_file = pd.read_csv(seed_path, sep="\t")

    if fasta_path is not None:
        fasta_file = ''
        fasta_pre_only_file = ''
        open(fasta_path, 'w').close()
        open(fasta_pre_only_path, 'w').close()

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
            overlaps = int(row['overlaps'])

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
                        seq5p_id += '|' + seed_file[seed_file['Seed'] == seq5p_seed]["Family"].iloc[0]
                    except:
                        seq5p_id += '|' + seq5p_seed

                if seq3p != '-':
                    seq3p_seed = seq3p[1:8].upper()
                    try:
                        seq3p_id += '|' + seed_file[seed_file['Seed'] == seq3p_seed]["Family"].iloc[0]
                    except:
                        seq3p_id += '|' + seq3p_seed

            if fasta_path is not None:
                if (seq5p != '-') & (mature_seq == 5):
                    fasta_file += f'>{seq5p_id}\n{seq5p}\n'
                    fasta_pre_only_file += f'>{seq5p_id}\n{hairpin}\n'

                if (seq3p != '-') & (mature_seq == 3):
                    fasta_file += f'>{seq3p_id}\n{seq3p}\n'
                    fasta_pre_only_file += f'>{seq3p_id}\n{hairpin}\n'

                if len(fasta_file) > 100000:
                    with open(fasta_path, 'a+') as f:
                        f.write(fasta_file)
                    fasta_file = ''

                if len(fasta_pre_only_file) > 100000:
                    with open(fasta_pre_only_path, 'a+') as f:
                        f.write(fasta_pre_only_file)
                    fasta_pre_only_file = ''

            start = int(positions.split('..')[0]) + 1
            end = int(positions.split('..')[1])
            if mature_seq == 5:
                seed = seq5p_id.split('|')[5]
                gff_row = [[f'{name}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={seq_id};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{overlaps}']]
            if mature_seq == 3:
                seed = seq3p_id.split('|')[5]
                gff_row = [[f'{name}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={seq_id};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{overlaps}']]
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
        with open(fasta_pre_only_path, 'a+') as f:
            f.write(fasta_pre_only_file)

    gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')
    gff3_pre_only.to_csv(output_pre_only, index=False, header=False, mode="a", sep='\t')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    output = None
    fasta_path = None
    seed_path = None
    species = None
    args = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '-o':
            output = sys.argv[i + 1]
        elif arg == '-seed':
            seed_path = sys.argv[i + 1]
        elif arg == '-s':
            species = sys.argv[i + 1]
        elif arg == '--create-fasta':
            fasta_path = sys.argv[i + 1]
        elif arg == '--exclude-c':
            exclude_c = sys.argv[i + 1]
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

    if not output:
        raise ('Output path is required (-o <path>)')


    run(output, fasta_path, seed_path)
