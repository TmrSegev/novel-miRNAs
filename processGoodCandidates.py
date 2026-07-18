#!/usr/bin/python
"""Generate unique_candidates CSV from debugging output (all species)."""

import argparse
import os

import pandas as pd

from pipeline_config import (
    candidates_csv_filename,
    candidates_dir,
    candidates_label,
    dev_condition,
    get_species_config,
)


def cluster_has_support(group_rows, cfg, tool_name, min_libraries=2):
    """Hofstenia: require ≥2 replicates of the same developmental condition."""
    if cfg["support_mode"] == "dev_condition_replicates":
        dev_counts = pd.Series(
            [dev_condition(cfg, row["Library"], tool_name) for row in group_rows]
        ).value_counts()
        return any(dev_counts >= 2)
    libraries = {row["Library"] for row in group_rows}
    return len(libraries) >= min_libraries


def _select_highest(current_group, precursor_col, mature_col, star_col, read_count_col):
    """Keep the highest-read-count row; mark sequence identity and cluster size."""
    best_idx = max(range(len(current_group)), key=lambda i: current_group[i][read_count_col])
    highest = current_group[best_idx].copy()
    same_precursor = len({x[precursor_col] for x in current_group}) == 1
    same_mature_star = len({(x[mature_col], x[star_col]) for x in current_group}) == 1
    highest["all_same"] = "yes" if same_precursor and same_mature_star else "no"
    highest["overlaps"] = len(current_group)
    filtered = [current_group[i] for i in range(len(current_group)) if i != best_idx]
    return highest, filtered


def process_group(current_group, cfg, tool_name, precursor_col, mature_col, star_col, read_count_col):
    """Pick one representative per ±20 bp cluster; optionally require multi-library support."""
    kept = []
    filtered = []
    if not current_group:
        return kept, filtered

    # Nematodes (unique_locus): keep one per cluster, including single-library loci.
    if cfg["support_mode"] == "unique_locus":
        highest, filtered = _select_highest(
            current_group, precursor_col, mature_col, star_col, read_count_col
        )
        kept.append(highest)
        return kept, filtered

    # Hofstenia: require ≥2 condition replicates; skip unsupported clusters.
    if len(current_group) <= 1:
        return kept, filtered
    if cluster_has_support(current_group, cfg, tool_name):
        highest, _ = _select_highest(
            current_group, precursor_col, mature_col, star_col, read_count_col
        )
        kept.append(highest)
    else:
        filtered.extend(current_group)
    return kept, filtered


def run(cfg, tool_name):
    scripts_dir = cfg["scripts_dir"]
    output_dir = candidates_dir(cfg)
    species = cfg["species"]
    label = candidates_label(cfg)
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

    kept_candidates = []
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
                kept, filt = process_group(
                    current_group, cfg, tool_name, precursor_col, mature_col, star_col, read_count_col
                )
                kept_candidates.extend(kept)
                filter_out.extend(filt)
                current_group = [row]
                first_start = row["start"]
        kept, filt = process_group(
            current_group, cfg, tool_name, precursor_col, mature_col, star_col, read_count_col
        )
        kept_candidates.extend(kept)
        filter_out.extend(filt)

    kept_df = pd.DataFrame(kept_candidates)
    filter_out_df = pd.DataFrame(filter_out)
    if cfg["support_mode"] == "dev_condition_replicates" and not kept_df.empty:
        kept_df["dev_condition"] = kept_df["Library"].apply(
            lambda lib: dev_condition(cfg, lib, tool_name)
        )
    if not kept_df.empty:
        tail = [c for c in ("dev_condition", "all_same", "overlaps") if c in kept_df.columns]
        head = [c for c in kept_df.columns if c not in tail]
        kept_df = kept_df[head + tail]
    kept_df.to_csv(
        os.path.join(output_dir, candidates_csv_filename(cfg, tool_name)), sep="\t", index=False
    )
    filter_out_df.to_csv(
        os.path.join(output_dir, f"{tool_name}_filterout.csv"), sep="\t", index=False
    )
    print(f"Processing complete for {tool_name}. Wrote {len(kept_df)} {label} candidates.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process debugging CSV into unique_candidates."
    )
    parser.add_argument("--tool", choices=["sRNAbench", "miRDeep"], required=True)
    parser.add_argument("-s", "--species", required=True, help="Species name")
    parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
    parser.add_argument("--variant", help='Genome variant, e.g. "new_genome" for alternate assembly track')
    args = parser.parse_args()
    cfg = get_species_config(args.species, args.base_path, variant=args.variant)
    run(cfg, args.tool)
