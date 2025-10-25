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



def run(output, fasta_path=None, seed_path=None, good_candidates=False, new_genome=False):
    """
    This Function will create GFF3 file from the sRNAbench output.
    :param output: base filename for the GFF3 output file.
    :param fasta_path: base filename for fasta files (if provided, creates fasta files).
    :param seed_path: a path to the seed file.
    :param good_candidates: if True, use pre-filtered good candidates.
    :param new_genome: if True, use new genome folder structure (Hofstenia_newGenome_*).
    :return:
    """
    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)
    gff3_pre_only = pd.DataFrame(columns=gff3_columns)
    # gff3_pre_only = gff3_pre_only.astype({"start": int})
    # gff3_pre_only = gff3_pre_only.astype({"end": int})
    
    # Determine base path based on genome version
    if new_genome:
        base_path = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/Hofstenia_newGenome/"
        output_dir = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia_newGenome/scripts/"
    else:
        base_path = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/sRNAtoolboxDB/out/"
        output_dir = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/scripts/"
    
    # Prepend output directory to all output file paths
    output = output_dir + output
    
    # Extract species name from output filename for output_pre_only
    # e.g., "Hofstenia_sRNAbench.gff3" -> "Hofstenia"
    species = output.split('/')[-1].split('_sRNAbench')[0]
    output_pre_only = output_dir + "{}_sRNAbench_pre_only.gff3".format(species)
    
    if fasta_path is not None:
        fasta_path = output_dir + fasta_path
        fasta_prefix = fasta_path.split('.fasta')[0]
        fasta_pre_only_path = fasta_prefix + "_pre_only.fasta"
        fasta_star_path = fasta_prefix + "_star.fasta"
    else:
        fasta_pre_only_path = None
        fasta_star_path = None

    # Uniting all remaining files, and all removed no find files
    table = None
    removed_no_find = None
    
    # Sample suffixes for all folders (folder names remain the same)
    sample_suffixes = ["EC1", "EC2", "EC3", "GA1", "GA2", "GA3", "DI1", "DI2", "DI3", 
                       "PDi1", "PDi2", "PDi3", "PDii1", "PDii2", "PDii3", "PL1", "PL2", "PL3", 
                       "PH1", "PH2", "PH3", "HL1", "HL2", "HL3", "IST1", "IST2", "IST3", 
                       "AMP1", "AMP2", "AMP3", "SMA1", "SMA2", "SMA3"]
    
    folders = ["Hofstenia_" + suffix for suffix in sample_suffixes]
    for folder in folders:
        to_add = pd.read_csv(base_path + folder + "/sRNAbench_remaining.csv", sep='\t')
        to_add_no_find = pd.read_csv(base_path + folder + "/sRNAbench_removed_no_find.csv", sep='\t')

        to_add["Library"] = folder

        if table is None:
            table = to_add
        else:
            table = pd.concat([table, to_add], ignore_index=True)

        if removed_no_find is None:
            removed_no_find = to_add_no_find
        else:
            removed_no_find = pd.concat([removed_no_find, to_add_no_find], ignore_index=True)

    removed_no_find.to_csv(output_dir + "all_sRNAbench_removed_no_find.csv", sep='\t', index=False)


    # Filtering by coordinates
    table = table.sort_values(['seqName', 'start', 'end'])
    table.to_csv(output_dir + 'debugging_Hofstenia_sRNAbench.csv', sep='\t', index=False)
    print("SAVED CSV")

    table['overlaps'] = np.zeros(len(table))
    no_overlaps = pd.DataFrame(columns=table.columns)

    for index, row in table.iterrows():
        if index in table.index:
            if len(row['hairpinSeq']) < 20:
                print(row)
            table['distance'] = (row['end'] - table['start']) / (row['end'] - row['start'])
            overlaps = table[(table['distance'] >= 0.6) & (table['distance'] <= 1)].tail(-1)
            overlaps = overlaps[overlaps['seqName'] == row['seqName']]
            table.loc[index, 'overlaps'] = len(overlaps)
            if len(overlaps) == 0:
                no_overlaps = no_overlaps.append(row)
                table = table.drop(index)
            else:
                table = table.drop(overlaps.index)
    print(table['overlaps'].value_counts().sort_index(ascending=False))
    if 'distance' in table.columns:
        table = table.drop(["distance"], axis=1)
    no_overlaps.to_csv(output_dir + 'removed_sRNAbench_no_overlaps.csv', sep='\t')
    # table.to_csv(output_dir + 'sRNAbench_all_remaining_filtered.csv', sep='\t', index=False)

    if good_candidates:
        if new_genome:
            good_candidates_path = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia_newGenome/good_candidates/sRNAbench_goodCandidates.csv"
        else:
            good_candidates_path = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/good_candidates/sRNAbench_goodCandidates.csv"
        try:
            table = pd.read_csv(good_candidates_path)
            if table.empty:
                print(f"Warning: Good candidates file is empty: {good_candidates_path}")
            table.to_csv(output_dir + 'sRNAbench_all_remaining_filtered.csv', sep='\t', index=False)
        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            print(f"Warning: Could not read good candidates file: {good_candidates_path}")
            print(f"Error: {e}")
            print("Continuing with filtered table from coordinate overlap filtering...")

    # Check if table is empty after all filtering
    if table.empty:
        print("Warning: No miRNA candidates remaining after filtering. Creating empty output files.")
        # Create empty output files
        with open(output, 'w') as file:
            file.write(version)
        with open(output_pre_only, 'w') as file:
            file.write(version)
        gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')
        gff3_pre_only.to_csv(output_pre_only, index=False, header=False, mode="a", sep='\t')
        if fasta_path is not None:
            open(fasta_path, 'w').close()
            open(fasta_pre_only_path, 'w').close()
            open(fasta_star_path, 'w').close()
        print("Empty output files created successfully.")
        return

    if seed_path:
        seed_file = pd.read_csv(seed_path, encoding='latin-1')

    if fasta_path is not None:
        fasta_file = ''
        fasta_pre_only_file = ''
        fasta_star_file = ''
        open(fasta_path, 'w').close()
        open(fasta_pre_only_path, 'w').close()
        open(fasta_star_path, 'w').close()

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
        overlaps = int(row['overlaps'])

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
                fasta_star_file += f'>{name5p}\n{seq3p}\n'
                fasta_pre_only_file += f'>{name5p}\n{hairpin}\n'
            if not pd.isnull(seq3p) and mature_seq == 3:
                fasta_file += f'>{name3p}\n{seq3p}\n'
                fasta_star_file += f'>{name3p}\n{seq5p}\n'
                fasta_pre_only_file += f'>{name3p}\n{hairpin}\n'

            if len(fasta_file) > 100000:
                with open(fasta_path, 'a+') as f:
                    f.write(fasta_file)
                fasta_file = ''

            if len(fasta_pre_only_file) > 100000:
                with open(fasta_pre_only_path, 'a+') as f:
                    f.write(fasta_pre_only_file)
                fasta_pre_only_file = ''

            if len(fasta_star_file) > 100000:
                with open(fasta_star_path, 'a+') as f:
                    f.write(fasta_star_file)
                fasta_star_file = ''

        if mature_seq == 5:
            seed = name5p.split('|')[4]
            gff_row = [[f'{seqId}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={name};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{origin};{overlaps}']]
        if mature_seq == 3:
            seed = name3p.split('|')[4]
            gff_row = [[f'{seqId}', '.', 'pre_miRNA', str(start), str(end), '.', strand, '.', f'ID={name};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{origin};{overlaps}']]
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
        with open(fasta_pre_only_path, 'a+') as f:
            f.write(fasta_pre_only_file)
        with open(fasta_star_path, 'a+') as f:
            f.write(fasta_star_file)

    gff3.to_csv(output, index=False, header=False, mode="a", sep='\t')
    gff3_pre_only.to_csv(output_pre_only, index=False, header=False, mode="a", sep='\t')


if __name__ == '__main__':
    output = None
    fasta_path = None
    seed_path = None
    species = None
    good_candidates = False
    new_genome = False
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
        elif arg == '--goodcandidates':
            good_candidates = sys.argv[i + 1]
            if good_candidates == "True":
                good_candidates = True
            else:
                good_candidates = False
        elif arg == '--new-genome':
            new_genome = sys.argv[i + 1]
            if new_genome == "True":
                new_genome = True
            else:
                new_genome = False
        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -o <filename>: output GFF3 filename (required).\n'
                  f' -s <species>: species name (required, used for output_pre_only filename).\n'
                  f' -seed <path> : classify the reads by seed file, should be separated by tab with columns.\n'
                  f' --create-fasta <filename>: create fasta file from the gff3 table.\n'
                  f' --goodcandidates <True/False>: use pre-filtered good candidates.\n'
                  f' --new-genome <True/False>: use new genome folder structure (Hofstenia_newGenome_*).\n'
                  )

            sys.exit()

    if not output:
        raise ('Output filename is required (-o <filename>)')
    if not species:
        raise ('Species name is required (-s <species>)')
    run(output, fasta_path, seed_path, good_candidates, new_genome)
# See PyCharm help at https://www.jetbrains.com/help/pycharm/