#!/usr/bin/env python3
import os
import argparse

# Paths
# base_dir = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Elegans/miRge"
# mature_in = os.path.join(base_dir, "all_candidates_mature.fasta")
# star_in   = os.path.join(base_dir, "all_candidates_star.fasta")
# combined_out = os.path.join(base_dir, "combined_mature_star.fasta")

def add_suffix_to_headers(in_path, suffix, out_handle):
    """
    Read a FASTA from in_path, add suffix (e.g. ";mature" or ";star")
    to every header line (>header), and write to out_handle.
    """
    with open(in_path) as f:
        for line in f:
            if line.startswith(">"):
                # strip trailing newline, add suffix, then newline
                name = line[1:].strip()
                new_header = f">{name}{suffix}\n"
                out_handle.write(new_header)
            else:
                out_handle.write(line)

def main():
    # --- Argument Parsing ---
    # Set up a parser to read command-line arguments
    parser = argparse.ArgumentParser(
        description="Combines 'all_candidates_mature.fasta' and 'all_candidates_star.fasta' from a base directory into a single FASTA file. It adds ';mature' and ';star' suffixes to the respective sequence headers."
    )
    # Add a required positional argument for the base directory path
    parser.add_argument(
        "--base_path",
        required=True,
        help="The path to the base directory containing the input FASTA files (all_candidates_mature.fasta and all_candidates_star.fasta)."
    )
    # Parse the arguments provided by the user
    args = parser.parse_args()
    base_dir = args.base_path

    # --- Path Construction ---
    # Construct the full paths for input and output files
    mature_in = os.path.join(base_dir, "all_candidates_mature_p.fasta")
    star_in = os.path.join(base_dir, "all_candidates_star_p.fasta")
    combined_out = os.path.join(base_dir, "combined_mature_star.fasta")

    # Open combined output for writing (will overwrite if exists)
    with open(combined_out, "w") as out_f:
        # Process mature file
        add_suffix_to_headers(mature_in, "", out_f)
        # Process star file
        add_suffix_to_headers(star_in, "", out_f)
    print(f"Wrote combined FASTA to {combined_out}")

if __name__ == "__main__":
    main()
