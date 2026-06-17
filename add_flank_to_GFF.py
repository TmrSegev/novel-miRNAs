import argparse
import os
import sys

from pipeline_config import get_species_config

parser = argparse.ArgumentParser(description="Add flanking regions to GFF3 pre_miRNA features")
parser.add_argument("-s", "--species", required=True, help="Species name (Elegans, Macrosperma, Sulstoni, Hofstenia)")
parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
parser.add_argument("--new-genome", action="store_true", help="Use Hofstenia_newGenome folder (Hofstenia only)")
args = parser.parse_args()

if args.new_genome and args.species == "Hofstenia":
    base_dir = os.path.join(args.base_path or "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq",
                            "Hofstenia_newGenome/scripts/")
    species_prefix = "Hofstenia"
else:
    cfg = get_species_config(args.species, args.base_path)
    base_dir = cfg["scripts_dir"]
    species_prefix = cfg["species"]

algorithms = ["mirdeep", "sRNAbench"]
FLANK = 10

for algo in algorithms:
    INPUT = os.path.join(base_dir, f"{species_prefix}_{algo}.gff3")
    OUTPUT = os.path.join(base_dir, f"{species_prefix}_{algo}_flanked_pre.gff3")

    print(f"Processing {INPUT} -> {OUTPUT}", file=sys.stderr)

    if not os.path.exists(INPUT):
        print(f"Warning: Input file not found: {INPUT}", file=sys.stderr)
        continue

    with open(INPUT) as check_file:
        lines = check_file.readlines()
        if all(line.startswith("#") for line in lines):
            print(f"Warning: Input file is empty (only headers): {INPUT}", file=sys.stderr)
            with open(OUTPUT, "w") as fout:
                fout.writelines(lines)
            continue

    with open(INPUT) as fin, open(OUTPUT, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            feature_type = cols[2]
            if feature_type == "pre_miRNA":
                start = int(cols[3]) - FLANK
                end = int(cols[4]) + FLANK
                cols[3] = str(max(start, 1))
                cols[4] = str(end)

                attributes = cols[8].split(";")

                if algo == "sRNAbench" and len(attributes) == 7:
                    attributes[4] = f"seed={attributes[4]}"
                    attributes[5] = f"novel={attributes[5]}"
                    attributes[6] = f"overlaps={attributes[6]}"
                    cols[8] = "|".join(attributes)
                elif algo == "mirdeep" and len(attributes) == 6:
                    attributes[4] = f"seed={attributes[4]}"
                    attributes[5] = f"overlaps={attributes[5]}"
                    cols[8] = "|".join(attributes)

            fout.write("\t".join(cols) + "\n")

print("Done.", file=sys.stderr)
