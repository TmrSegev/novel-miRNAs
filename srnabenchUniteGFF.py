#!/usr/bin/python
"""Unite per-library sRNAbench results and write GFF3 (all species)."""

import argparse
import os

import pandas as pd

from pipeline_config import get_species_config, load_seed_table, resolve_seed_path, srnabench_folder_name
from pipeline_overlap import deduplicate_coordinate_overlaps

ids_dic = {}


def handle_given_name(name, df, column):
    if len(df[df[column] == name]) > 1:
        if name not in ids_dic:
            ids_dic[name] = 0
        ids_dic[name] += 1
        name = f"{name}_{ids_dic[name]}"
    return name


def run(output, species_name, fasta_path=None, seed_path=None, good_candidates=False,
        debug_only=False, base_path=None, variant=None):
    if debug_only and good_candidates:
        raise ValueError("Cannot use --debug-only together with --goodcandidates True")

    cfg = get_species_config(species_name, base_path, variant=variant)
    output_dir = cfg["scripts_dir"]
    srnabench_out = cfg["srnabench_out_dir"]
    good_candidates_dir = cfg["good_candidates_dir"]
    libraries = cfg["libraries"]
    species = cfg["display_species"]

    version = "##gff-version 3\n"
    gff3_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    output = os.path.join(output_dir, output)
    output_pre_only = os.path.join(output_dir, f"{species}_sRNAbench_pre_only.gff3")

    if fasta_path is not None:
        fasta_path = os.path.join(output_dir, fasta_path)
        fasta_prefix = fasta_path.split(".fasta")[0]
        fasta_pre_only_path = fasta_prefix + "_pre_only.fasta"
        fasta_star_path = fasta_prefix + "_star.fasta"
    else:
        fasta_pre_only_path = None
        fasta_star_path = None

    good_candidates_path = os.path.join(good_candidates_dir, "sRNAbench_goodCandidates.csv")
    if good_candidates and not debug_only and os.path.exists(good_candidates_path):
        table = pd.read_csv(good_candidates_path, sep="\t")
        if table.empty:
            print(f"Warning: Good candidates file is empty: {good_candidates_path}")
        table.to_csv(os.path.join(output_dir, "sRNAbench_all_remaining_filtered.csv"), sep="\t", index=False)
        print(f"Loaded good candidates from {good_candidates_path} (skipping per-library re-unite)")
        _write_srnabench_outputs(
            table, output, output_pre_only, fasta_path, fasta_pre_only_path, fasta_star_path,
            version, gff3_columns, cfg, seed_path,
        )
        return

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
        for path in (output, output_pre_only):
            with open(path, "w") as fh:
                fh.write(version)
        return

    removed_no_find.to_csv(os.path.join(output_dir, "all_sRNAbench_removed_no_find.csv"), sep="\t", index=False)

    table = table.sort_values(["seqName", "start", "end"])
    debug_path = os.path.join(output_dir, f"debugging_{species}_sRNAbench.csv")
    table.to_csv(debug_path, sep="\t", index=False)
    print(f"Saved debugging CSV: {debug_path}")

    if debug_only:
        print("Debug-only mode: stopping after debugging CSV (run processGoodCandidates.py next).")
        return

    table, no_overlaps = deduplicate_coordinate_overlaps(table, "seqName")
    print(table["overlaps"].value_counts().sort_index(ascending=False))
    no_overlaps.to_csv(os.path.join(output_dir, "removed_sRNAbench_no_overlaps.csv"), sep="\t")

    if good_candidates:
        try:
            table = pd.read_csv(good_candidates_path, sep="\t")
            if table.empty:
                print(f"Warning: Good candidates file is empty: {good_candidates_path}")
            table.to_csv(os.path.join(output_dir, "sRNAbench_all_remaining_filtered.csv"), sep="\t", index=False)
        except (FileNotFoundError, pd.errors.EmptyDataError) as e:
            print(f"Warning: Could not read good candidates file: {good_candidates_path}")
            print(f"Error: {e}")

    _write_srnabench_outputs(
        table, output, output_pre_only, fasta_path, fasta_pre_only_path, fasta_star_path,
        version, gff3_columns, cfg, seed_path,
    )


def _write_srnabench_outputs(table, output, output_pre_only, fasta_path, fasta_pre_only_path,
                             fasta_star_path, version, gff3_columns, cfg, seed_path):
    gff3 = pd.DataFrame(columns=gff3_columns)
    gff3_pre_only = pd.DataFrame(columns=gff3_columns)

    if table.empty:
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
    for _, row in table.iterrows():
        intersection_index += 1
        name = handle_given_name(row["name"], table, "name")
        seq_id = row["seqName"]
        name5p = handle_given_name(row["5pname"], table, "5pname")
        seq5p = row["5pseq"]
        name3p = handle_given_name(row["3pname"], table, "3pname")
        seq3p = row["3pseq"]
        strand = row["strand"]
        hairpin = row["hairpinSeq"]
        start = row["start"]
        end = row["end"]
        origin = row["origin"]
        overlaps = int(row["overlaps"])

        if row["5pRC"] >= row["3pRC"]:
            name5p += "|m"
            name3p += "|s"
            mature_arm = 5
            rc_mature = row["5pRC"]
            rc_star = row["3pRC"]
        else:
            name5p += "|s"
            name3p += "|m"
            mature_arm = 3
            rc_mature = row["3pRC"]
            rc_star = row["5pRC"]

        seq5p_freq = len(table[(table["5pseq"] == seq5p) | (table["3pseq"] == seq5p)])
        seq3p_freq = len(table[(table["5pseq"] == seq3p) | (table["3pseq"] == seq3p)])
        name5p += f"|{seq5p_freq}"
        name3p += f"|{seq3p_freq}"
        name5p += f"|index={intersection_index}"
        name3p += f"|index={intersection_index}"

        if seed_file is not None:
            if not pd.isnull(seq5p):
                seq5p_seed = seq5p[1:8].upper().replace("T", "U")
                try:
                    name5p += "|" + seed_file[seed_file["Seed"] == seq5p_seed]["Family"].iloc[0]
                except IndexError:
                    name5p += "|" + seq5p_seed
            if not pd.isnull(seq3p):
                seq3p_seed = seq3p[1:8].upper().replace("T", "U")
                try:
                    name3p += "|" + seed_file[seed_file["Seed"] == seq3p_seed]["Family"].iloc[0]
                except IndexError:
                    name3p += "|" + seq3p_seed

        if fasta_path is not None:
            if not pd.isnull(seq5p) and mature_arm == 5:
                fasta_file += f">{name5p}\n{seq5p}\n"
                fasta_star_file += f">{name5p}\n{seq3p}\n"
                fasta_pre_only_file += f">{name5p}\n{hairpin}\n"
            if not pd.isnull(seq3p) and mature_arm == 3:
                fasta_file += f">{name3p}\n{seq3p}\n"
                fasta_star_file += f">{name3p}\n{seq5p}\n"
                fasta_pre_only_file += f">{name3p}\n{hairpin}\n"

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

        if mature_arm == 5:
            seed = name5p.split("|")[4]
        else:
            seed = name3p.split("|")[4]
        gff_row = [[seq_id, ".", "pre_miRNA", str(start), str(end), ".", strand, ".",
                    f"ID={name};RC_m={rc_mature};RC_s={rc_star};index={intersection_index};{seed};{origin};{overlaps}"]]
        gff3_pre_only = pd.concat([gff3_pre_only, pd.DataFrame(gff_row, columns=gff3_columns)], ignore_index=True)

        if strand == "+":
            try:
                offset5p = len(hairpin.split(seq5p)[0])
                gff_row.append([seq_id, ".", "miRNA", start + offset5p, start + offset5p + len(seq5p) - 1,
                                ".", strand, ".", f"ID={name5p}"])
            except Exception:
                pass
            try:
                offset3p = len(hairpin.split(seq3p)[0])
                gff_row.append([seq_id, ".", "miRNA", start + offset3p, start + offset3p + len(seq3p) - 1,
                                ".", strand, ".", f"ID={name3p}"])
            except Exception:
                pass
        else:
            try:
                offset5p = len(hairpin.split(seq5p)[0])
                gff_row.append([seq_id, ".", "miRNA", end - offset5p - len(seq5p) + 1, end - offset5p,
                                ".", strand, ".", f"ID={name5p}"])
            except Exception:
                pass
            try:
                offset3p = len(hairpin.split(seq3p)[0])
                gff_row.append([seq_id, ".", "miRNA", end - offset3p - len(seq3p) + 1, end - offset3p,
                                ".", strand, ".", f"ID={name3p}"])
            except Exception:
                pass

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
    parser = argparse.ArgumentParser(description="Unite per-library sRNAbench results and write GFF3.")
    parser.add_argument("-o", required=True, help="Output GFF3 filename")
    parser.add_argument("-s", required=True, dest="species", help="Species name")
    parser.add_argument("-seed", help="Seed family file (default from species config)")
    parser.add_argument("--create-fasta", dest="fasta_path", help="Output mature FASTA filename")
    parser.add_argument("--goodcandidates", default="False", help="True to use good_candidates CSV")
    parser.add_argument("--debug-only", action="store_true",
                        help="Write debugging CSV only; skip GFF/FASTA (Step A of good_candidates workflow)")
    parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
    parser.add_argument("--variant", help='Genome variant, e.g. "new_genome" for alternate assembly track')
    args = parser.parse_args()
    run(
        args.o,
        args.species,
        fasta_path=args.fasta_path,
        seed_path=args.seed,
        good_candidates=args.goodcandidates == "True",
        debug_only=args.debug_only,
        base_path=args.base_path,
        variant=args.variant,
    )
