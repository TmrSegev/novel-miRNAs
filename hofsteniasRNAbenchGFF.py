#!/usr/bin/python

import sys
import pandas as pd
import numpy as np

ids_dic = {}


def handleGivenName(name, df, column):
    """
    This Function handle the given name of miRNA,
    name of miRNA used as a key, most be unique.
    :param name: sRNAbench given name.
    :param df: dataframe which including all the records.
    :param column: column name.
    :return: unique name.
    """
    if len(df[df[column] == name]) > 1:
        if name not in ids_dic:
            ids_dic[name] = 0
        ids_dic[name] += 1
        name = f'{name}_{ids_dic[name]}'
    return name



def run(output, fasta_path=None, seed_path=None):
    """
    This Function will create GFF3 file from the sRNAbench output.
    :param seed_path: a path to the seed file.
    :param fasta_path: a path to create fasta file from the gff3 table.
    :param additional: additonal sRNAbench output prediction file.
    :param input: sRNAbench output prediction files 'novel.txt', 'novel454.txt'
    :param output: output path of the GFF formatted file.
    :return:
    """
    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)
    gff3_pre_only = pd.DataFrame(columns=gff3_columns)
    # gff3_pre_only = gff3_pre_only.astype({"start": int})
    # gff3_pre_only = gff3_pre_only.astype({"end": int})
    output_pre_only = "{}_sRNAbench_pre_only.gff3".format(species)

    # Uniting all remaining files
    table = None
    folders = ["Hofstenia_EC1", "Hofstenia_EC2", "Hofstenia_EC3", "Hofstenia_GA1", "Hofstenia_GA2", "Hofstenia_GA3", "Hofstenia_DI1", "Hofstenia_DI2", "Hofstenia_DI3", "Hofstenia_PDi1", "Hofstenia_PDi2", "Hofstenia_PDi3", "Hofstenia_PDii1", "Hofstenia_PDii2", "Hofstenia_PDii3", "Hofstenia_PL1", "Hofstenia_PL2", "Hofstenia_PL3", "Hofstenia_PH1", "Hofstenia_PH2", "Hofstenia_PH3", "Hofstenia_HL1", "Hofstenia_HL2", "Hofstenia_HL3", "Hofstenia_IST1", "Hofstenia_IST2", "Hofstenia_IST3", "Hofstenia_AMP1", "Hofstenia_AMP2", "Hofstenia_AMP3", "Hofstenia_SMA1", "Hofstenia_SMA2", "Hofstenia_SMA3"]
    for folder in folders:
        to_add = pd.read_csv("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/" + folder + "/sRNAbench_remaining.csv", sep='\t')
        if table is None:
            table = to_add
        else:
            table = pd.concat([table, to_add], ignore_index=True)

    # Filtering by coordinates
    table = table.sort_values(['seqName', 'start', 'end'])

    for index, row in table.iterrows():
        table['overlaps'] = (row['end'] - table['start'])/(row['end'] - row['start'])
        overlaps = table[(table['overlaps'] >= 0.6) & (table['overlaps'] <= 1)].tail(-1)
        overlaps = overlaps[overlaps['seqName'] == row['seqName']]
        table = table.drop(overlaps.index)


    if seed_path:
        seed_file = pd.read_csv(seed_path, sep='\t')

    if fasta_path is not None:
        fasta_file = ''
        open(fasta_path, 'w').close()

    intersection_index = -1 # Used later to intersect the table with miRdeep, blast and featurecounts results.
    for index, row in table.iterrows():
        intersection_index += 1
        name = handleGivenName(row['name'], table, 'name')
        seqId = row['seqName']
        name5p = handleGivenName(row['5pname'], table, '5pname')
        seq5p = row['5pseq']
        name3p = handleGivenName(row['3pname'], table, '3pname')
        seq3p = row['3pseq']
        strand = row['strand']
        hairpin = row['hairpinSeq']
        start = row['start']
        end = row['end']
        origin = row['origin']

        if row['5pRC'] >= row['3pRC']:
            name5p += '|m'
            name3p += '|s'
            mature_seq = 5
            rc_mature = row['5pRC']
            rc_star = row['3pRC']
        else:
            name5p += '|s'
            name3p += '|m'
            mature_seq = 3
            rc_mature = row['3pRC']
            rc_star = row['5pRC']

        seq5p_freq = len(table[(table['5pseq'] == seq5p) | (table['3pseq'] == seq5p)])
        seq3p_freq = len(table[(table['5pseq'] == seq3p) | (table['3pseq'] == seq3p)])

        name5p += f'|{seq5p_freq}'
        name3p += f'|{seq3p_freq}'

        name5p += f'|index={intersection_index}'
        name3p += f'|index={intersection_index}'

        if seed_path is not None:
            if not pd.isnull(seq5p):
                seq5p_seed = seq5p[1:8].upper().replace("T", "U")
                try:
                    name5p += '|' + seed_file[seed_file['Seed'] == seq5p_seed]["Family"].iloc[0]
                except:
                    name5p += '|' + seq5p_seed

            if not pd.isnull(seq3p):
                seq3p_seed = seq3p[1:8].upper().replace("T", "U")
                try:
                    name3p += '|' + seed_file[seed_file['Seed'] == seq3p_seed]["Family"].iloc[0]
                except:
                    name3p += '|' + seq3p_seed

        if fasta_path is not None:
            if not pd.isnull(seq5p) and mature_seq == 5:
                fasta_file += f'>{name5p}\n{seq5p}\n'
            if not pd.isnull(seq3p) and mature_seq == 3:
                fasta_file += f'>{name3p}\n{seq3p}\n'

            if len(fasta_file) > 100000:
                with open(fasta_path, 'a+') as f:
                    f.write(fasta_file)
                fasta_file = ''

        if mature_seq == 5:
            seed = name5p.split('|')[4]
            gff_row = [[f'{seqId}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={name};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{origin}']]
        if mature_seq == 3:
            seed = name3p.split('|')[4]
            gff_row = [[f'{seqId}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={name};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{origin}']]
        gff3_pre_only = gff3_pre_only.append(gff_row)

        if strand == '+':
            try:
                offset5p = len(hairpin.split(seq5p)[0])
                start5p = start + offset5p
                end5p = start + offset5p + len(seq5p) - 1
                gff_row.append([seqId, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={name5p}'])
            except:
                pass

            try:
                offset3p = len(hairpin.split(seq3p)[0])
                start3p = start + offset3p
                end3p = start + offset3p + len(seq3p) - 1
                gff_row.append([seqId, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={name3p}'])
            except:
                pass

        else:
            try:
                offset5p = len(hairpin.split(seq5p)[0])
                end5p = end - offset5p
                start5p = end - offset5p - len(seq5p) + 1
                gff_row.append([seqId, '.', 'miRNA', start5p, end5p, '.', strand, '.', f'ID={name5p}'])
            except:
                pass

            try:
                offset3p = len(hairpin.split(seq3p)[0])
                end3p = end - offset3p
                start3p = end - offset3p - len(seq3p) + 1
                gff_row.append([seqId, '.', 'miRNA', start3p, end3p, '.', strand, '.', f'ID={name3p}'])
            except:
                pass

        miRNAs = pd.DataFrame(gff_row, columns=gff3_columns)

        gff3 = gff3.append(miRNAs)

    with open(output, 'w') as file:
        file.write(version)

    with open(output_pre_only, 'w') as file:
        file.write(version)

    if fasta_path is not None:
        with open(fasta_path, 'a+') as f:
            f.write(fasta_file)

    gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')
    gff3_pre_only.to_csv(output_pre_only, index=False, header=False, mode="a", sep='\t')


if __name__ == '__main__':
    output = None
    fasta_path = None
    seed_path = None
    species = None
    args = []
    for i in range(1, len(sys.argv), 2):
        arg = sys.argv[i]
        if arg == '-o':
            output = sys.argv[i + 1]
        elif arg == '-seed':
            seed_path = sys.argv[i + 1]
        elif arg == '--create-fasta':
            fasta_path = sys.argv[i + 1]
        elif arg == '-s':
            species = sys.argv[i + 1]
        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i <path>: sRNAbench prediction output, like novel.txt/novel451.txt.\n'
                  f' -a <path>: additional input file.\n'
                  f' -o <path>: output path.\n'
                  f' -seed <path> : classify the reads by seed file, should be separated by tab with columns.\n'
                  f' --create-fasta <path>: create fasta file from the gff3 table.\n'
                  )

            sys.exit()

    if not output:
        raise ('Output path is required (-o <path>)')
    run(output, fasta_path, seed_path)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/