#!/usr/bin/python

import sys
import pandas as pd
import numpy as np


def filterNovel451(novel451, novel):
    """
    :param novel451: dataframe of the novel451.txt additional input.
    :param novel: dataframe of the novel.txt input.
    """
    deleted_input = pd.DataFrame(columns=novel451.columns)
    for index_451, row_451 in novel451.iterrows():
        seq3p = row_451['3pseq']
        seq5p = row_451['5pseq']
        if max(row_451['3pRC'], row_451['5pRC']) < threshold:
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

def filterNovel(novel):
    """
    :param novel451: dataframe of the novel451.txt additional input.
    :param novel: dataframe of the novel.txt input.
    """
    deleted_input = pd.DataFrame(columns=novel.columns)
    for index, row in novel.iterrows():
        if max(row['3pRC'], row['5pRC']) < threshold:
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


def run(input, additional=None):
    """
    This Function will create GFF3 file from the sRNAbench output.
    :param seed_path: a path to the seed file.
    :param fasta_path: a path to create fasta file from the gff3 table.
    :param additional: additonal sRNAbench output prediction file.
    :param input: sRNAbench output prediction files 'novel.txt', 'novel454.txt'
    :param output: output path of the GFF formatted file.
    :return:
    """
    table = pd.read_csv(input, sep='\t')
    table["origin"] = "novel"
    table, deleted_input = filterNovel(table)

    if additional:
        table_to_add = pd.read_csv(additional, sep='\t')
        table_to_add["origin"] = "novel451"
        table_to_add, table_to_delete = filterNovel451(table_to_add, table)
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

    table.to_csv('sRNAbench_remaining.csv', sep='\t', index=False)
    deleted_input.to_csv('sRNAbench_removed.csv', sep='\t', index = False)

if __name__ == '__main__':
    input = None
    add = None
    species = None
    threshold = 10
    args = []
    for i in range(1, len(sys.argv), 2):
        arg = sys.argv[i]
        if arg == '-i':
            input = sys.argv[i + 1]
        elif arg == '-a':
            add = sys.argv[i + 1]
        elif arg == '-s':
            species = sys.argv[i + 1]
        elif arg == '-t':
            threshold = sys.argv[i + 1]
        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -i <path>: sRNAbench prediction output, like novel.txt/novel451.txt.\n'
                  f' -a <path>: additional input file.\n'
                  f' -o <path>: output path.\n'
                  f' -seed <path> : classify the reads by seed file, should be separated by tab with columns.\n'
                  f' --create-fasta <path>: create fasta file from the gff3 table.\n'
                  )

            sys.exit()

    if not input:
        raise ('Input path is required (-i <path>)')
    run(input, add)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/