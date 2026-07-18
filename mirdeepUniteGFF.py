#!/usr/bin/python
"""Unite per-library miRDeep results and write GFF3 (all species)."""

import argparse
import os

import pandas as pd

from pipeline_config import (
    candidates_csv_path,
    candidates_label,
    get_species_config,
    load_seed_table,
    mirdeep_folder_name,
    resolve_seed_path,
)
from pipeline_overlap import deduplicate_coordinate_overlaps


def get_seq_id(row):
    try:
        return row["provisional id"]
    except KeyError:
        seq_id = row["mature miRBase miRNA"]
        return seq_id.replace("-3p", "").replace("-5p", "")


def run(output, species_name, fasta_path=None, seed_path=None, unique_candidates=False,
        debug_only=False, base_path=None, variant=None):
    if debug_only and unique_candidates:
        raise ValueError("Cannot use --debug-only together with --uniquecandidates True")

    cfg = get_species_config(species_name, base_path, variant=variant)
    output_dir = cfg["scripts_dir"]
    mirdeep_out = cfg["mirdeep_out_dir"]
    libraries = cfg["libraries"]
    species = cfg["display_species"]
    cand_label = candidates_label(cfg)
    curated_path = candidates_csv_path(cfg, "miRDeep")

    version = "##gff-version 3\n"
    gff3_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    output = os.path.join(output_dir, output)
    output_pre_only = os.path.join(output_dir, f"{species}_mirdeep_pre_only.gff3")

    if fasta_path is not None:
        fasta_path = os.path.join(output_dir, fasta_path)
        fasta_prefix = fasta_path.split(".fasta")[0]
        fasta_pre_only_path = fasta_prefix + "_pre_only.fasta"
        fasta_star_path = fasta_prefix + "_star.fasta"
    else:
        fasta_pre_only_path = None
        fasta_star_path = None

    if unique_candidates and not debug_only and os.path.exists(curated_path):
        unique_df = pd.read_csv(curated_path, sep="\t")
        if unique_df.empty:
            print(f"Warning: {cand_label.capitalize()} candidates file is empty: {curated_path}")
        unique_df.to_csv(os.path.join(output_dir, "mirdeep_all_remaining_filtered.csv"), sep="\t", index=False)
        print(f"Loaded {cand_label} candidates from {curated_path} (skipping per-library re-unite)")
        _write_mirdeep_outputs(
            [unique_df], output, output_pre_only, fasta_path, fasta_pre_only_path, fasta_star_path,
            version, gff3_columns, cfg, seed_path,
        )
        return

    filtered_input = []
    for file_num in range(1, 3):
        table = None
        for library in libraries:
            folder = mirdeep_folder_name(cfg, library)
            remaining_path = os.path.join(mirdeep_out, folder, f"remaining_file_{file_num}.csv")
            if not os.path.exists(remaining_path):
                print(f"Warning: missing {remaining_path}, skipping")
                continue
            to_add = pd.read_csv(remaining_path, sep="\t")
            to_add["Library"] = folder
            table = to_add if table is None else pd.concat([table, to_add], ignore_index=True)

        if table is None or table.empty:
            print(f"Table for iteration {file_num} is empty - skipping filtering.")
            pd.DataFrame().to_csv(
                os.path.join(output_dir, f"removed_mirdeep_{file_num}_no_overlaps.csv"), sep="\t"
            )
            continue

        table = table.sort_values(["precursor coordinate"])
        debug_path = os.path.join(output_dir, f"debugging_{species}_miRDeep_{file_num}.csv")
        table.to_csv(debug_path, sep="\t", index=False)
        print(f"Saved debugging CSV: {debug_path}")

        if debug_only:
            continue

        table["chr"] = table["precursor coordinate"].str.split(":", expand=True)[0]
        table["positions"] = table["precursor coordinate"].str.split(":", expand=True)[1].astype(str)
        table["start"] = table["positions"].str.split(".", expand=True)[0].astype(int)
        table["end"] = table["positions"].str.split(".", expand=True)[2].astype(int)

        table, no_overlaps = deduplicate_coordinate_overlaps(table, "chr")
        print(table["overlaps"].value_counts().sort_index(ascending=False))
        filtered_input.append(table)
        no_overlaps.to_csv(os.path.join(output_dir, f"removed_mirdeep_{file_num}_no_overlaps.csv"), sep="\t")
        table = table.rename(
            {
                "tag id": "provisional id",
                "estimated probability that the miRNA is a true positive":
                    "estimated probability that the miRNA candidate is a true positive",
            },
            axis=1,
        )

    if debug_only:
        print("Debug-only mode: stopping after debugging CSV (run processGoodCandidates.py next).")
        return

    if unique_candidates:
        try:
            unique_df = pd.read_csv(curated_path, sep="\t")
            if unique_df.empty:
                print(f"Warning: {cand_label.capitalize()} candidates file is empty: {curated_path}")
            unique_df.to_csv(os.path.join(output_dir, "mirdeep_all_remaining_filtered.csv"), sep="\t", index=False)
            filtered_input = [unique_df]
        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            print(f"Warning: Could not read {cand_label} candidates file: {curated_path}")
            print(f"Error: {e}")
            print("Continuing with filtered inputs from coordinate overlap filtering...")

    _write_mirdeep_outputs(
        filtered_input, output, output_pre_only, fasta_path, fasta_pre_only_path, fasta_star_path,
        version, gff3_columns, cfg, seed_path,
    )


def _write_mirdeep_outputs(filtered_input, output, output_pre_only, fasta_path, fasta_pre_only_path,
                           fasta_star_path, version, gff3_columns, cfg, seed_path):
    gff3 = pd.DataFrame(columns=gff3_columns)
    gff3_pre_only = pd.DataFrame(columns=gff3_columns)

    if not filtered_input or all(df.empty for df in filtered_input):
        print("Warning: No miRNA candidates remaining after filtering. Creating empty output files.")
        for path in (output, output_pre_only):
            with open(path, "w") as fh:
                fh.write(version)
        gff3.to_csv(output, index=False, header=False, mode="a", sep="\t")
        gff3_pre_only.to_csv(output_pre_only, index=False, header=False, mode="a", sep="\t")
        if fasta_path is not None:
            for p in (fasta_path, fasta_pre_only_path, fasta_star_path):
                open(p, "w").close()
        return

    seed_file = None
    resolved_seed = resolve_seed_path(cfg, seed_path)
    if resolved_seed and os.path.exists(resolved_seed):
        seed_file = load_seed_table(resolved_seed, cfg["seed_encoding"], cfg.get("seed_sep", "\t"))

    if fasta_path is not None:
        fasta_file = ""
        fasta_pre_only_file = ""
        fasta_star_file = ""
        open(fasta_path, "w").close()
        open(fasta_pre_only_path, "w").close()
        open(fasta_star_path, "w").close()

    intersection_index = -1
    for input_df in filtered_input:
        for _, row in input_df.iterrows():
            intersection_index += 1
            details = row["precursor coordinate"]
            name = details.split(":")[0]
            positions = details.split(":")[1]
            strand = details.split(":")[2]
            seq_id = get_seq_id(row)
            star_seq = row["consensus star sequence"]
            mature_seq = row["consensus mature sequence"]
            hairpin = row["consensus precursor sequence"]
            rc_mature = row["mature read count"]
            rc_star = row["star read count"]
            overlaps = int(row["overlaps"])

            star_position = hairpin.index(star_seq)
            mature_position = hairpin.index(mature_seq)

            seq5p_id = seq_id + "|5p"
            seq3p_id = seq_id + "|3p"

            if star_position > mature_position:
                seq5p = row["consensus mature sequence"]
                seq3p = row["consensus star sequence"]
                seq5p_freq = len(input_df[input_df["consensus mature sequence"] == seq5p])
                seq3p_freq = len(input_df[input_df["consensus star sequence"] == seq3p])
                seq5p_id += f"|m|{seq5p_freq}"
                seq3p_id += f"|s|{seq3p_freq}"
                mature_arm = 5
            else:
                seq5p = row["consensus star sequence"]
                seq3p = row["consensus mature sequence"]
                seq5p_freq = len(input_df[input_df["consensus star sequence"] == seq5p])
                seq3p_freq = len(input_df[input_df["consensus mature sequence"] == seq3p])
                seq5p_id += f"|s|{seq5p_freq}"
                seq3p_id += f"|m|{seq3p_freq}"
                mature_arm = 3

            seq5p_id += f"|index={intersection_index}"
            seq3p_id += f"|index={intersection_index}"

            if seed_file is not None:
                if seq5p != "-":
                    seq5p_seed = seq5p[1:8].upper()
                    try:
                        seq5p_id += "|" + seed_file[seed_file["Seed"] == seq5p_seed]["Family"].iloc[0]
                    except IndexError:
                        seq5p_id += "|" + seq5p_seed
                if seq3p != "-":
                    seq3p_seed = seq3p[1:8].upper()
                    try:
                        seq3p_id += "|" + seed_file[seed_file["Seed"] == seq3p_seed]["Family"].iloc[0]
                    except IndexError:
                        seq3p_id += "|" + seq3p_seed

            if fasta_path is not None:
                if (seq5p != "-") & (mature_arm == 5):
                    fasta_file += f">{seq5p_id}\n{seq5p}\n"
                    fasta_pre_only_file += f">{seq5p_id}\n{hairpin}\n"
                    fasta_star_file += f">{seq5p_id}\n{seq3p}\n"
                if (seq3p != "-") & (mature_arm == 3):
                    fasta_file += f">{seq3p_id}\n{seq3p}\n"
                    fasta_pre_only_file += f">{seq3p_id}\n{hairpin}\n"
                    fasta_star_file += f">{seq3p_id}\n{seq5p}\n"

                if len(fasta_file) > 100000:
                    with open(fasta_path, "a+") as fh:
                        fh.write(fasta_file)
                    fasta_file = ""
                if len(fasta_pre_only_file) > 100000:
                    with open(fasta_pre_only_path, "a+") as fh:
                        fh.write(fasta_pre_only_file)
                    fasta_pre_only_file = ""
                if len(fasta_star_file) > 100000:
                    with open(fasta_star_path, "a+") as fh:
                        fh.write(fasta_star_file)
                    fasta_star_file = ""

            start = int(positions.split("..")[0]) + 1
            end = int(positions.split("..")[1])
            if mature_arm == 5:
                seed = seq5p_id.split("|")[5]
            else:
                seed = seq3p_id.split("|")[5]
            gff_row = [[name, ".", "pre_miRNA", str(start), str(end), ".", strand, ".",
                        f"ID={seq_id};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{overlaps}"]]
            gff3_pre_only = pd.concat([gff3_pre_only, pd.DataFrame(gff_row, columns=gff3_columns)], ignore_index=True)

            if strand == "+":
                if seq5p != "-":
                    offset5p = len(hairpin.split(seq5p)[0])
                    gff_row.append([name, ".", "miRNA", start + offset5p, start + offset5p + len(seq5p) - 1,
                                    ".", strand, ".", f"ID={seq5p_id}"])
                if seq3p != "-":
                    offset3p = len(hairpin.split(seq3p)[0])
                    gff_row.append([name, ".", "miRNA", start + offset3p, start + offset3p + len(seq3p) - 1,
                                    ".", strand, ".", f"ID={seq3p_id}"])
            else:
                if seq5p != "-":
                    offset5p = len(hairpin.split(seq5p)[0])
                    gff_row.append([name, ".", "miRNA", end - offset5p - len(seq5p) + 1, end - offset5p,
                                    ".", strand, ".", f"ID={seq5p_id}"])
                if seq3p != "-":
                    offset3p = len(hairpin.split(seq3p)[0])
                    gff_row.append([name, ".", "miRNA", end - offset3p - len(seq3p) + 1, end - offset3p,
                                    ".", strand, ".", f"ID={seq3p_id}"])

            gff3 = pd.concat([gff3, pd.DataFrame(gff_row, columns=gff3_columns)], ignore_index=True)

    with open(output, "w") as fh:
        fh.write(version)
    with open(output_pre_only, "w") as fh:
        fh.write(version)

    if fasta_path is not None:
        with open(fasta_path, "a+") as fh:
            fh.write(fasta_file)
        with open(fasta_pre_only_path, "a+") as fh:
            fh.write(fasta_pre_only_file)
        with open(fasta_star_path, "a+") as fh:
            fh.write(fasta_star_file)

    gff3.to_csv(output, index=False, header=False, mode="a", sep="\t")
    gff3_pre_only.to_csv(output_pre_only, index=False, header=False, mode="a", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unite per-library miRDeep results and write GFF3.")
    parser.add_argument("-o", required=True, help="Output GFF3 filename")
    parser.add_argument("-s", required=True, dest="species", help="Species name")
    parser.add_argument("-seed", help="Seed family file (default from species config)")
    parser.add_argument("--create-fasta", dest="fasta_path", help="Output mature FASTA filename")
    parser.add_argument(
        "--uniquecandidates",
        default="False",
        help="True to use unique_candidates CSV",
    )
    parser.add_argument(
        "--debug-only",
        action="store_true",
        help="Write debugging CSV only; skip GFF/FASTA (Step A of unique_candidates workflow)",
    )
    parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
    parser.add_argument("--variant", help='Genome variant, e.g. "new_genome" for alternate assembly track')
    args = parser.parse_args()
    run(
        args.o,
        args.species,
        fasta_path=args.fasta_path,
        seed_path=args.seed,
        unique_candidates=args.uniquecandidates == "True",
        debug_only=args.debug_only,
        base_path=args.base_path,
        variant=args.variant,
    )
