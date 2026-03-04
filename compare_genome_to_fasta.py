import gffutils
from Bio import SeqIO
import os
import pandas as pd
import argparse

# --- Argument Parsing ---
# Set up a parser to read command-line arguments
parser = argparse.ArgumentParser(
    description=""
)
# Add a required positional argument for the base directory path
parser.add_argument("--species", help="")
parser.add_argument("--dir", help="")
parser.add_argument("--genome_fasta", help="")
parser.add_argument("--gff", help="")
parser.add_argument("--premir", help="")
parser.add_argument("--mature_star", help="")
parser.add_argument("--output", help="")
parser.add_argument("--db", help="")
# Parse the arguments provided by the user
args = parser.parse_args()
SPECIES = args.species
dir = args.dir
genome_fasta = args.genome_fasta
gff_file = args.gff
pre_miRNA_fasta = args.premir
mature_star = args.mature_star
output = args.output
db = args.db

# File paths
# oscar_dir = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge/"
# genome_fasta = "/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/Hmia.030120.fasta"
# gff_file = oscar_dir + "miRNA_250203_reformatted.gff3"
# pre_miRNA_fasta = oscar_dir + "pre_miR_1050_no_pre_in_seqid.fa"
# mature_star_fasta = oscar_dir + "combined_mature_star_1050.fa"
# output_file = oscar_dir + "genome_to_fasta_comparison_new.csv"
# db_file = oscar_dir + "miRNA_250203_reformatted.db"

gff_file = os.path.join(dir, gff_file)
pre_miRNA_fasta = os.path.join(dir, pre_miRNA_fasta)
mature_star_fasta = os.path.join(dir, mature_star)
output_file = os.path.join(dir, output)
db_file = os.path.join(dir, output)

# Load genome sequence
genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# Load reference sequences
def load_reference_sequences(fasta_file):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

pre_miRNA_ref = load_reference_sequences(pre_miRNA_fasta)
mature_star_ref = load_reference_sequences(mature_star_fasta)

# Initialize result list
results = []

# Function to extract sequence from genome
def extract_sequence(seq_id, start, end, strand):
    sequence = genome[seq_id].seq[start-1:end]
    return str(sequence.reverse_complement()) if strand == '-' else str(sequence)

# Compare sequences and save results
def compare_sequences(feature_type, extracted_seq, ref_seq, feature_id, strand):
    match_status = "Match" if extracted_seq == ref_seq else "Mismatch"
    results.append([feature_id, feature_type, extracted_seq, ref_seq, strand, match_status])

# Check if the database file exists; if not, create it
if not os.path.exists(db_file):
    print("Creating GFF database...")
    db = gffutils.create_db(gff_file, db_file, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)
else:
    db = gffutils.FeatureDB(db_file)

# Parse GFF and compare sequences
for feature in db.all_features():
    if feature.featuretype in ["miRNA_primary_transcript", "miRNA"]:
        feature_id = feature.attributes.get("ID", ["Unknown"])[0]
        start, end, strand = feature.start, feature.end, feature.strand
        seq_id = feature.seqid

        extracted_seq = extract_sequence(seq_id, start, end, strand)

        # Compare with corresponding reference sequence
        if feature.featuretype == "miRNA_primary_transcript" and feature_id in pre_miRNA_ref:
            compare_sequences("miRNA_primary_transcript", extracted_seq, pre_miRNA_ref[feature_id], feature_id, strand)
        elif feature.featuretype == "miRNA" and feature_id in mature_star_ref:
            compare_sequences("miRNA", extracted_seq, mature_star_ref[feature_id], feature_id, strand)

# Convert results to DataFrame and save as CSV
df = pd.DataFrame(results, columns=["Feature ID", "Feature Type", "Extracted Sequence", "Reference Sequence", "Strand", "Match Status"])
df.to_csv(output_file, index=False)

print("Processing complete. Results saved in:", output_file)

