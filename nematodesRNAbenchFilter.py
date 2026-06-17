#!/usr/bin/python
"""Per-library sRNAbench filter for nematode (and configurable) pipelines."""

import argparse
import os
import re
import sys

import numpy as np
import pandas as pd
from Bio import pairwise2

from pipeline_config import DEFAULT_NCRNA_DIR, get_species_config

threshold = 10


def filterNovel451(novel451):
    novel451 = novel451.copy()
    novel451["Removal Reason"] = "novel451"
    deleted_input = novel451.copy()
    return pd.DataFrame(columns=novel451.columns), deleted_input


def filterNovel(novel):
    deleted_input = pd.DataFrame(columns=novel.columns)
    for index, row in novel.iterrows():
        if max(row["3pRC"], row["5pRC"]) < threshold:
            row = row.copy()
            row["Removal Reason"] = "Weak mature signal"
            deleted_input = deleted_input.append(row)
            novel.drop(index=index, inplace=True)
            continue
        if row["matureBindings"] < 14:
            row = row.copy()
            row["Removal Reason"] = "Hairpin does not have enough pairings"
            deleted_input = deleted_input.append(row)
            novel.drop(index=index, inplace=True)
    return novel, deleted_input


def start_5p(row):
    if row["5pseq"] != "nan":
        for a in pairwise2.align.globalms(
            row["hairpinSeq"], row["5pseq"], 1, -1, -1, -1, penalize_end_gaps=False
        ):
            return re.search(r"[ATCG]", a.seqB).start()
    return 0


def end_3p(row):
    if row["3pseq"] != "nan":
        for a in pairwise2.align.globalms(
            row["hairpinSeq"], row["3pseq"], 1, -1, -1, -1, penalize_end_gaps=False
        ):
            return re.search(r"[ATCG]", a.seqB).start() + len(row["3pseq"])
    return len(row["hairpinSeq"])


def cut_hairpin(row):
    return row["hairpinSeq"][row["start_5p"] : row["end_3p"]]


def filter_ncrna(table, deleted_input, ncrna_dir):
    table = table.copy()
    table["5pseq"] = table["5pseq"].astype(str)
    table["3pseq"] = table["3pseq"].astype(str)
    nc_files = {
        "rRNA": "Caenorhabditis_rRNA.fasta",
        "snoRNA": "Caenorhabditis_snoRNA.fasta",
        "snRNA": "Caenorhabditis_snRNA.fasta",
        "tRNA": "Caenorhabditis_tRNA.fasta",
    }
    for index, row in table.iterrows():
        for label, fname in nc_files.items():
            path = os.path.join(ncrna_dir, fname)
            if not os.path.exists(path):
                continue
            with open(path) as f:
                content = f.read()
            if (row["5pseq"] in content) or (row["3pseq"] in content):
                row = row.copy()
                row["Removal Reason"] = label
                deleted_input = deleted_input.append(row)
                table.drop(index=index, inplace=True)
                break
    return table, deleted_input


def run(input_path, additional=None, ncrna_dir=DEFAULT_NCRNA_DIR):
    table = pd.read_csv(input_path, sep="\t")
    table["origin"] = "novel"

    if table.empty:
        pd.DataFrame(columns=["Removal Reason"]).to_csv("sRNAbench_removed.csv", sep="\t", index=False)
        table.to_csv("sRNAbench_remaining.csv", sep="\t", index=False)
        pd.DataFrame(columns=table.columns).to_csv("sRNAbench_removed_no_find.csv", sep="\t", index=False)
        return

    table, deleted_input = filterNovel(table)

    if additional and os.path.exists(additional):
        table_to_add = pd.read_csv(additional, sep="\t")
        table_to_add["origin"] = "novel451"
        table_to_add, table_to_delete = filterNovel451(table_to_add)
        deleted_input = deleted_input.append(table_to_delete)
        table = table.append(table_to_add)
    elif additional:
        print(f"Warning: Additional file not found: {additional}. Skipping novel451.")

    table, deleted_input = filter_ncrna(table, deleted_input, ncrna_dir)

    table["start_5p"] = table.apply(start_5p, axis=1)
    table["end_3p"] = table.apply(end_3p, axis=1)
    remove_no_find = table[(table["start_5p"] == -1) | (table["end_3p"] == -1)]
    remove_no_find.to_csv("sRNAbench_removed_no_find.csv", sep="\t", index=False)
    table = table[(table["start_5p"] != -1) & (table["end_3p"] != -1)]
    table["hairpinSeq"] = table.apply(cut_hairpin, axis=1)

    table.loc[table["strand"] == "+", "end"] = table["start"] + table["end_3p"] - 1
    table.loc[table["strand"] == "+", "start"] = table["start"] + table["start_5p"]
    table.loc[table["strand"] == "-", "start"] = table["end"] - table["end_3p"] + 1
    table.loc[table["strand"] == "-", "end"] = table["end"] - table["start_5p"]
    table = table.drop(["start_5p", "end_3p"], axis=1)

    table.to_csv("sRNAbench_remaining.csv", sep="\t", index=False)
    deleted_input.to_csv("sRNAbench_removed.csv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter per-library sRNAbench output.")
    parser.add_argument("-i", required=True, help="novel.txt path")
    parser.add_argument("-a", help="novel451.txt path")
    parser.add_argument(
        "-t", "--filter-mc", type=int, default=10,
        help="Filter if max(5pRC,3pRC) < threshold (nematodes: use 100)",
    )
    parser.add_argument("--ncrna-dir", default=DEFAULT_NCRNA_DIR)
    args = parser.parse_args()
    threshold = args.filter_mc
    run(args.i, args.a, args.ncrna_dir)
