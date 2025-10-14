import sys

# A list of the algorithms whose output will be processed.
algorithms = ["mirdeep", "sRNAbench"]
# The number of bases to add to each side of the pre_miRNA feature.
FLANK = 10

for algo in algorithms:
    # These are placeholder paths.
    # You should replace them with the actual paths to your files.
    INPUT = f"/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_{algo}.gff3"
    OUTPUT = f"/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/scripts/Hofstenia_{algo}_flanked_pre.gff3"

    print(f"Processing {INPUT} -> {OUTPUT}", file=sys.stderr)

    with open(INPUT) as fin, open(OUTPUT, "w") as fout:
        for line in fin:
            # Preserve comment/header lines verbatim.
            if line.startswith("#"):
                fout.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            feature_type = cols[2]
            print(feature_type)
            if feature_type == "pre_miRNA":
                # --- 1. Existing logic to add flanks (unchanged) ---
                start = int(cols[3]) - FLANK
                end = int(cols[4]) + FLANK
                # GFF positions are 1-based, so start cannot be less than 1.
                cols[3] = str(max(start, 1))
                cols[4] = str(end)

                # --- 2. MODIFIED LOGIC to fix the malformed attribute column ---
                attributes = cols[8].split(';')

                # Handle the sRNAbench format, which has 7 attribute fields.
                if algo == "sRNAbench" and len(attributes) == 7:
                    attributes[4] = f"seed={attributes[4]}"
                    attributes[5] = f"novel={attributes[5]}"
                    attributes[6] = f"overlaps={attributes[6]}"

                    # Join the list back into a valid GFF3 attribute string.
                    cols[8] = "|".join(attributes)

                # Handle the mirdeep format, which has 6 attribute fields.
                elif algo == "mirdeep" and len(attributes) == 6:
                    # The 5th field is a known miRNA ID or a novel mature sequence.
                    # 'name' is a general and descriptive key for this value.
                    attributes[4] = f"seed={attributes[4]}"
                    # The 6th field is a numeric value, which we'll label as 'score'.
                    attributes[5] = f"overlaps={attributes[5]}"

                    # Join the list back into a valid GFF3 attribute string.
                    cols[8] = "|".join(attributes)

            # Write the modified or original line to the output file.
            fout.write("\t".join(cols) + "\n")

print("Done.", file=sys.stderr)