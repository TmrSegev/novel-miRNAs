#!/usr/bin/python
"""Generate good_candidates CSV from debugging output (all species)."""

import argparse
import os

import pandas as pd

from pipeline_config import dev_condition, get_species_config


def cluster_has_support(group_rows, cfg, tool_name, min_libraries=2):
    if cfg["support_mode"] == "dev_condition_replicates":
        dev_counts = pd.Series(
            [dev_condition(cfg, row["Library"], tool_name) for row in group_rows]
        ).value_counts()
        return any(dev_counts >= 2)
    libraries = {row["Library"] for row in group_rows}
    return len(libraries) >= min_libraries


def process_group(current_group, cfg, tool_name, precursor_col, mature_col, star_col, read_count_col):
    good = []
    filtered = []
    if len(current_group) <= 1:
        return good, filtered
    if cluster_has_support(current_group, cfg, tool_name):
        highest = max(current_group, key=lambda x: x[read_count_col])
        same_precursor = len({x[precursor_col] for x in current_group}) == 1
        same_mature_star = len({(x[mature_col], x[star_col]) for x in current_group}) == 1
        highest = highest.copy()
        highest["all_same"] = "yes" if same_precursor and same_mature_star else "no"
        highest["overlaps"] = len(current_group)
        good.append(highest)
    else:
        filtered.extend(current_group)
    return good, filtered


def run(cfg, tool_name):
    scripts_dir = cfg["scripts_dir"]
    output_dir = cfg["good_candidates_dir"]
    species = cfg["species"]
    os.makedirs(output_dir, exist_ok=True)

    if tool_name == "miRDeep":
        file1 = os.path.join(scripts_dir, f"debugging_{species}_miRDeep_1.csv")
        file2 = os.path.join(scripts_dir, f"debugging_{species}_miRDeep_2.csv")
        df1 = pd.read_csv(file1, sep="\t")
        df2 = pd.read_csv(file2, sep="\t", names=df1.columns, skiprows=1)
        df = pd.concat([df1, df2], ignore_index=True)
        combined_debug = os.path.join(scripts_dir, f"debugging_{species}_miRDeep.csv")
        df.to_csv(combined_debug, sep="\t", index=False)
    else:
        input_file = os.path.join(scripts_dir, f"debugging_{species}_sRNAbench.csv")
        df = pd.read_csv(input_file, sep="\t")
        df = df[~df["origin"].str.contains("novel451", na=False)]

    dev_condition_col = "Library"
    precursor_col = "hairpinSeq" if tool_name == "sRNAbench" else "consensus precursor sequence"
    mature_col = "3pseq" if tool_name == "sRNAbench" else "consensus mature sequence"
    star_col = "5pseq" if tool_name == "sRNAbench" else "consensus star sequence"
    read_count_col = "totalRC" if tool_name == "sRNAbench" else "total read count"

    if tool_name == "miRDeep":
        df[["scaffold", "coordinates", "strand"]] = df["precursor coordinate"].str.split(":", expand=True)
        df[["start", "end"]] = df["coordinates"].str.split("\\.\\.", expand=True)
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
        df = df.sort_values(by=["scaffold", "strand", "start"]).copy()
        df.drop(columns=["coordinates"], inplace=True)
    else:
        df["scaffold"] = df["seqName"]
        df["strand"] = df["strand"]
        df = df.sort_values(by=["scaffold", "strand", "start"]).copy()

    good_candidates = []
    filter_out = []

    for (scaffold, strand), sub_df in df.groupby(["scaffold", "strand"]):
        current_group = []
        first_start = None
        for _, row in sub_df.iterrows():
            row = row.copy()
            if not current_group:
                first_start = row["start"]
                current_group.append(row)
            elif row["start"] <= first_start + 20:
                current_group.append(row)
            else:
                good, filt = process_group(
                    current_group, cfg, tool_name, precursor_col, mature_col, star_col, read_count_col
                )
                good_candidates.extend(good)
                filter_out.extend(filt)
                current_group = [row]
                first_start = row["start"]
        good, filt = process_group(
            current_group, cfg, tool_name, precursor_col, mature_col, star_col, read_count_col
        )
        good_candidates.extend(good)
        filter_out.extend(filt)

    good_candidates_df = pd.DataFrame(good_candidates)
    filter_out_df = pd.DataFrame(filter_out)
    if cfg["support_mode"] == "dev_condition_replicates" and not good_candidates_df.empty:
        good_candidates_df["dev_condition"] = good_candidates_df["Library"].apply(
            lambda lib: dev_condition(cfg, lib, tool_name)
        )
    if not good_candidates_df.empty:
        tail = [c for c in ("dev_condition", "all_same", "overlaps") if c in good_candidates_df.columns]
        head = [c for c in good_candidates_df.columns if c not in tail]
        good_candidates_df = good_candidates_df[head + tail]
    good_candidates_df.to_csv(
        os.path.join(output_dir, f"{tool_name}_goodCandidates.csv"), sep="\t", index=False
    )
    filter_out_df.to_csv(
        os.path.join(output_dir, f"{tool_name}_filterout.csv"), sep="\t", index=False
    )
    print(f"Processing complete for {tool_name}. Wrote {len(good_candidates_df)} good candidates.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process debugging CSV into good_candidates.")
    parser.add_argument("--tool", choices=["sRNAbench", "miRDeep"], required=True)
    parser.add_argument("-s", "--species", required=True, help="Species name")
    parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
    parser.add_argument("--variant", help='Genome variant, e.g. "new_genome" for alternate assembly track')
    args = parser.parse_args()
    cfg = get_species_config(args.species, args.base_path, variant=args.variant)
    run(cfg, args.tool)
