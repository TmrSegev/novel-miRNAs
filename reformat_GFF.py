import pandas as pd
import argparse

# --- Argument Parsing ---
# Set up a parser to read command-line arguments
parser = argparse.ArgumentParser(
    description=""
)
# Add a required positional argument for the base directory path
parser.add_argument("--input", help="")
parser.add_argument("--output", help="")
parser.add_argument("--oscar",  action='store_true', help="")
# Parse the arguments provided by the user
args = parser.parse_args()

input_file = args.input
output_file = args.output
oscar = args.oscar
debug_counter = 0
if not oscar:
    print("OSCAR FALSE")


# Define input and output file paths
# input_file = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge/miRNA.250203.gff"  # Update with your actual file path
# output_file = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge/miRNA_250203_reformatted.gff3"

# Define GFF column headers
gff_columns = [
    "seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"
]

# Load the file into a DataFrame
miRNA_df = pd.read_csv(
    input_file,
    sep="\t",
    header=None,
    names=gff_columns,
    comment="#",  # Skip comment lines
    dtype=str  # Ensure all data is treated as strings to preserve formatting
)

# Replace "pre_miRNA" with "miRNA_primary_transcript"
miRNA_df["feature"] = miRNA_df["feature"].replace("pre_miRNA", "miRNA_primary_transcript")

# Define a function to reformat attributes and remove "-pre" suffixes
def reformat_attributes(row):
    attributes = row["attributes"].replace("-pre", "")
    print("attributes:")
    print(attributes)
    feature = row["feature"]
    attributes_list = attributes.split(";")
    print("attributes list:")
    print(attributes_list)

    gene_id = None
    if oscar:
        gene_id = next((x.split("=")[1] for x in attributes_list if x.startswith("gene_id")), None)
    else:
        print("Gene id:")
        print(attributes_list[0].strip().startswith("ID="))
        gene_id = next((x.split("=", 1)[1] for x in attributes_list if x.strip().startswith("ID=")), None)

    if feature == "miRNA_primary_transcript":
        return f"ID={gene_id};Alias={gene_id};Name={gene_id}"
    elif feature == "miRNA":
        parent = None
        if oscar:
            parent = next((x.split("=")[1] for x in attributes_list if x.startswith("Parent=")), None)
        else:
            # CORRECTED: Added .strip() to handle leading whitespace
            parent = next((x.split("=", 1)[1] for x in attributes_list if x.strip().startswith("Derives_from=")), None)
        return f"ID={gene_id};Alias={gene_id};Name={gene_id};Derives_from={parent}"

    return attributes

# Apply reformatting
miRNA_df["attributes"] = miRNA_df.apply(reformat_attributes, axis=1)

# Write to output GFF3 file
with open(output_file, "w") as f:
    miRNA_df.to_csv(f, sep="\t", header=False, index=False, line_terminator="\n")  # Write data

print(f"Reformatted file saved as {output_file}")
