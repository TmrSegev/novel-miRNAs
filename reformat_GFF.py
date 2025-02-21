import pandas as pd

# Define input and output file paths
input_file = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/miRge/miRNA.250203.gff"  # Update with your actual file path
output_file = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/miRge/miRNA_250203_reformatted.gff3"

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
    attributes = row["attributes"].replace("-pre", "")  # Remove "-pre" suffixes
    feature = row["feature"]
    gene_id = next((x.split("=")[1] for x in attributes.split(";") if x.startswith("gene_id")), None)

    if feature == "miRNA_primary_transcript":
        return f"ID={gene_id};Alias={gene_id};Name={gene_id}"
    elif feature == "miRNA":
        parent = next((x.split("=")[1] for x in attributes.split(";") if x.startswith("Parent")), None)
        return f"ID={gene_id};Alias={gene_id};Name={gene_id};Derives_from={parent}"

    return attributes

# Apply reformatting
miRNA_df["attributes"] = miRNA_df.apply(reformat_attributes, axis=1)

# Write to output GFF3 file
with open(output_file, "w") as f:
    miRNA_df.to_csv(f, sep="\t", header=False, index=False, line_terminator="\n")  # Write data

print(f"Reformatted file saved as {output_file}")
