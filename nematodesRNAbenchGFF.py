#!/usr/bin/python

import argparse
import os
import sys

import numpy as np
import pandas as pd

from pipeline_config import get_species_config, srnabench_folder_name
from pipeline_overlap import deduplicate_coordinate_overlaps

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



def run(output, species_name, fasta_path=None, seed_path=None, good_candidates=False, base_path=None):
    cfg = get_species_config(species_name, base_path)
    output_dir = cfg["scripts_dir"]
    srnabench_out = cfg["srnabench_out_dir"]
    good_candidates_dir = cfg["good_candidates_dir"]
    libraries = cfg["libraries"]
    species = cfg["species"]

    version = "##gff-version 3\n"
    gff3_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff3 = pd.DataFrame(columns=gff3_columns)
    gff3_pre_only = pd.DataFrame(columns=gff3_columns)

    output = os.path.join(output_dir, output)
    output_pre_only = os.path.join(output_dir, f"{species}_sRNAbench_pre_only.gff3")
    
    if fasta_path is not None:
        fasta_path = os.path.join(output_dir, fasta_path)
        fasta_prefix = fasta_path.split('.fasta')[0]
        fasta_pre_only_path = fasta_prefix + "_pre_only.fasta"
        fasta_star_path = fasta_prefix + "_star.fasta"
    else:
        fasta_pre_only_path = None
        fasta_star_path = None

    table = None
    removed_no_find = None

    for library in libraries:
        folder = srnabench_folder_name(cfg, library)
        remaining_path = os.path.join(srnabench_out, folder, "sRNAbench_remaining.csv")
        no_find_path = os.path.join(srnabench_out, folder, "sRNAbench_removed_no_find.csv")
        if not os.path.exists(remaining_path):
            print(f"Warning: missing {remaining_path}, skipping library {library}")
            continue
        to_add = pd.read_csv(remaining_path, sep="\t")
        to_add_no_find = pd.read_csv(no_find_path, sep="\t") if os.path.exists(no_find_path) else pd.DataFrame()
        to_add["Library"] = folder
        table = to_add if table is None else pd.concat([table, to_add], ignore_index=True)
        if removed_no_find is None:
            removed_no_find = to_add_no_find
        elif not to_add_no_find.empty:
            removed_no_find = pd.concat([removed_no_find, to_add_no_find], ignore_index=True)

    if table is None or table.empty:
        print("Warning: No sRNAbench remaining files found. Creating empty outputs.")
        with open(output, "w") as file:
            file.write(version)
        with open(output_pre_only, "w") as file:
            file.write(version)
        return

    removed_no_find.to_csv(os.path.join(output_dir, "all_sRNAbench_removed_no_find.csv"), sep="\t", index=False)

    table = table.sort_values(["seqName", "start", "end"])
    debug_path = os.path.join(output_dir, f"debugging_{species}_sRNAbench.csv")
    table.to_csv(debug_path, sep="\t", index=False)
    print(f"Saved debugging CSV: {debug_path}")

    table, no_overlaps = deduplicate_coordinate_overlaps(table, "seqName")
    print(table["overlaps"].value_counts().sort_index(ascending=False))
    no_overlaps.to_csv(os.path.join(output_dir, "removed_sRNAbench_no_overlaps.csv"), sep="\t")

    if good_candidates:
        good_candidates_path = os.path.join(good_candidates_dir, "sRNAbench_goodCandidates.csv")
        try:
            table = pd.read_csv(good_candidates_path, sep="\t")
            if table.empty:
                print(f"Warning: Good candidates file is empty: {good_candidates_path}")
            table.to_csv(os.path.join(output_dir, "sRNAbench_all_remaining_filtered.csv"), sep="\t", index=False)
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
        gff3_pre_only = pd.concat([gff3_pre_only, pd.DataFrame(gff_row, columns=gff3_columns)], ignore_index=True)

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

        gff3 = pd.concat([gff3, miRNAs], ignore_index=True)

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unite per-library sRNAbench results and write GFF3.")
    parser.add_argument("-o", required=True, help="Output GFF3 filename")
    parser.add_argument("-s", required=True, dest="species", help="Species name (Elegans, Macrosperma, Sulstoni)")
    parser.add_argument("-seed", help="Seed family CSV path")
    parser.add_argument("--create-fasta", dest="fasta_path", help="Output mature FASTA filename")
    parser.add_argument("--goodcandidates", default="False", help="True to use good_candidates CSV")
    parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
    args = parser.parse_args()
    run(
        args.o,
        args.species,
        fasta_path=args.fasta_path,
        seed_path=args.seed,
        good_candidates=args.goodcandidates == "True",
        base_path=args.base_path,
    )