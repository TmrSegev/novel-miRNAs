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


def filterNovel451(novel451, novel, threshold_mc):
    """
    :param novel451: dataframe of the novel451.txt additional input.
    :param novel: dataframe of the novel.txt input.
    """
    deleted_input = pd.DataFrame(columns=novel451.columns)
    for index_451, row_451 in novel451.iterrows():
        seq3p = row_451['3pseq']
        seq5p = row_451['5pseq']
        if max(row_451['3pRC'], row_451['5pRC']) < threshold_mc:
            row_451["Removal Reason"] = 'Weak mature signal'
            deleted_input = deleted_input.append(row_451)
            novel451.drop(index_451, inplace=True)
            continue
        if row_451['matureBindings'] < 14:
            row_451["Removal Reason"] = 'Hairpin does not have enough pairings'
            deleted_input = deleted_input.append(row_451)
            novel451.drop(index_451, inplace=True)
            continue
        # This segment removes mature seq from novel 451 if it is also present in novel.
        for index_novel, row_novel in novel.iterrows():
            if((row_novel["3pseq"] == seq5p) | (row_novel["5pseq"] == seq3p) | (row_novel["3pseq"] == seq3p) | (row_novel["5pseq"] == seq5p)):
                row_451["Removal Reason"] = 'Duplicate in novel, originally novel451'
                deleted_input = deleted_input.append(row_451)
                novel451.drop(index_451, inplace=True)
                break
    return novel451, deleted_input

def filterNovel(novel, threshold_mc):
    """
    :param novel451: dataframe of the novel451.txt additional input.
    :param novel: dataframe of the novel.txt input.
    """
    deleted_input = pd.DataFrame(columns=novel.columns)
    for index, row in novel.iterrows():
        if max(row['3pRC'], row['5pRC']) < threshold_mc:
            row["Removal Reason"] = 'Weak mature signal'
            deleted_input = deleted_input.append(row)
            novel.drop(index, inplace=True)
            continue
        if row['matureBindings'] < 14:
            row["Removal Reason"] = 'Hairpin does not have enough pairings'
            deleted_input = deleted_input.append(row)
            novel.drop(index, inplace=True)
            continue
    return novel, deleted_input


def start_5p(row):
    if row['5pseq'] != "nan":
        return row['hairpinSeq'].find(row['5pseq'])
    else:
        return 0


def end_3p(row):
    if row['3pseq'] != "nan":
        return row['hairpinSeq'].find(row['3pseq']) + len(row['3pseq'])
    else:
        return len(row['hairpinSeq'])


def cut_hairpin(row):
    return row['hairpinSeq'][row['start_5p']:row['end_3p']]


def run(input, output, threshold_mc, additional=None, fasta_path=None, seed_path=None):
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
    table = pd.read_csv(input, sep='\t')
    table["origin"] = "novel"
    table, deleted_input = filterNovel(table, threshold_mc)

    if seed_path:
        seed_file = pd.read_csv(seed_path, sep='\t')

    if fasta_path is not None:
        fasta_file = ''
        open(fasta_path, 'w').close()

    if additional:
        table_to_add = pd.read_csv(additional, sep='\t')
        table_to_add["origin"] = "novel451"
        table_to_add, table_to_delete = filterNovel451(table_to_add, table, threshold_mc)
        deleted_input = deleted_input.append(table_to_delete)
        table = table.append(table_to_add)

    # Filter non coding RNA
    table['5pseq'] = table['5pseq'].astype(str)
    table['3pseq'] = table['3pseq'].astype(str)
    for index, row in table.iterrows():
        with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_rRNA.fasta') as f:
            if (row['5pseq'] in f.read()) or (row['3pseq'] in f.read()):
                row['Removal Reason'] = 'rRNA'
                deleted_input = deleted_input.append(row)
                table.drop(index=index, inplace=True)
                f.close()
        with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_snoRNA.fasta') as f:
            if (row['5pseq'] in f.read()) or (row['3pseq'] in f.read()):
                row['Removal Reason'] = 'snoRNA'
                deleted_input = deleted_input.append(row)
                table.drop(index=index, inplace=True)
                f.close()
        with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_snRNA.fasta') as f:
            if (row['5pseq'] in f.read()) or (row['3pseq'] in f.read()):
                row['Removal Reason'] = 'snRNA'
                deleted_input = deleted_input.append(row)
                table.drop(index=index, inplace=True)
                f.close()
        with open('/sise/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/ncRNAs_Caenorhabditis/Caenorhabditis_tRNA.fasta') as f:
            if (row['5pseq'] in f.read()) or (row['3pseq'] in f.read()):
                row['Removal Reason'] = 'tRNA'
                deleted_input = deleted_input.append(row)
                table.drop(index=index, inplace=True)
                f.close()

    # Trim the sequences and adjust start/end
    # table['5pseq'] = table['5pseq'].fillna(np.nan).replace([np.nan], [None])
    # table['3pseq'] = table['3pseq'].fillna(np.nan).replace([np.nan], [None])

    table['start_5p'] = table.apply(lambda row : start_5p(row), axis=1)

    table['end_3p'] = table.apply(lambda row : end_3p(row), axis=1)

    table['hairpinSeq'] = table.apply(lambda row : cut_hairpin(row), axis=1)
    # table['hairpinSeq'] = table['hairpinSeq'].str[table['start_5p']:table['end_3p']]
    table['end'] = table['start'] + table['end_3p'] - 1
    table['start'] = table['start'] + table['start_5p']
    table = table.drop(['start_5p', 'end_3p'], axis=1)

    table.to_csv('sRNAbench_remaining.csv', sep='\t')
    deleted_input.to_csv('sRNAbench_removed.csv', sep='\t')

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
    input = None
    output = None
    add = None
    fasta_path = None
    seed_path = None
    species = None
    threshold_mc = None
    args = []
    for i in range(1, len(sys.argv), 2):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        elif arg == '-o':
            output = sys.argv[i + 1]
        elif arg == '-a':
            add = sys.argv[i + 1]
        elif arg == '-seed':
            seed_path = sys.argv[i + 1]
        elif arg == '--create-fasta':
            fasta_path = sys.argv[i + 1]
        elif arg == '-s':
            species = sys.argv[i + 1]
        elif arg == '--filter-mc':
            threshold_mc = sys.argv[i + 1]
        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i <path>: sRNAbench prediction output, like novel.txt/novel451.txt.\n'
                  f' -a <path>: additional input file.\n'
                  f' -o <path>: output path.\n'
                  f' -seed <path> : classify the reads by seed file, should be separated by tab with columns.\n'
                  f' --create-fasta <path>: create fasta file from the gff3 table.\n'
                  f' --filter-mc <int>: Filter if max(mature read count, star read count) < threshold of filter-mc.\n'
                  )

            sys.exit()

    if not input:
        raise ('Input path is required (-i <path>)')
    if not output:
        raise ('Output path is required (-o <path>)')

    if threshold_mc is not None:
        threshold_mc = float(threshold_mc)

    run(input, output, threshold_mc, add, fasta_path, seed_path)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/