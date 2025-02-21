# import gffutils
# from Bio import SeqIO
# import os
#
# # File paths
# oscar_dir = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/miRge/"
# genome_fasta = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/Hmia.030120.fasta"
# gff_file = oscar_dir + "miRNA_250203_reformatted.gff3"
# # gff_file = oscar_dir + "miRNA_no_pre_reformatted_v4.gff3"
# pre_miRNA_fasta = oscar_dir + "pre_miR_1050_no_pre_in_seqid.fa"
# mature_star_fasta = oscar_dir + "combined_mature_star_1050.fa"
# output_file = oscar_dir + "genome_to_fasta_comparison_new.txt"
# db_file = oscar_dir + "miRNA_250203_reformatted.db"
# # db_file = oscar_dir + "miRNA_reformatted_v4.db"
#
# # Load genome sequence
# genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
#
# # Load reference sequences
# def load_reference_sequences(fasta_file):
#     return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
#
# pre_miRNA_ref = load_reference_sequences(pre_miRNA_fasta)
# mature_star_ref = load_reference_sequences(mature_star_fasta)
#
# # Initialize counters
# match_count = 0
# mismatch_count = 0
#
# # Function to extract sequence from genome
# def extract_sequence(seq_id, start, end, strand):
#     sequence = genome[seq_id].seq[start-1:end]
#     return sequence.reverse_complement() if strand == '-' else sequence
#
# # Compare sequences and save mismatches to file
# def compare_sequences(feature_type, extracted_seq, ref_seq, feature_id, output_handle, strand):
#     global match_count, mismatch_count
#     if str(extracted_seq) != ref_seq:
#         mismatch_count += 1
#         output_handle.write(f"Mismatch in {feature_type} ({feature_id}):\n")
#         output_handle.write(f"Genome: {extracted_seq}\n")
#         output_handle.write(f"Reference: {ref_seq}\n")
#         output_handle.write(f"Strand: {strand}\n\n")
#     else:
#         match_count += 1
#
# # Check if the database file exists; if not, create it
# if not os.path.exists(db_file):
#     print("Creating GFF database...")
#     db = gffutils.create_db(gff_file, db_file, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)
# else:
#     db = gffutils.FeatureDB(db_file)
#
# # Parse GFF and compare sequences
# with open(output_file, "w") as output_handle:
#     for feature in db.all_features():
#         if feature.featuretype in ["miRNA_primary_transcript", "miRNA"]:
#
#             feature_id = feature.attributes.get("ID", ["Unknown"])[0]
#             start, end, strand = feature.start, feature.end, feature.strand
#             seq_id = feature.seqid
#
#             print(f"Feature ID: {feature_id}, Feature Type: {feature.featuretype}")
#             print(f"Available pre_miRNA_ref keys (first 5): {list(pre_miRNA_ref.keys())[:5]}")
#             print(f"Available mature_star_ref keys (first 5): {list(mature_star_ref.keys())[:5]}")
#
#             extracted_seq = extract_sequence(seq_id, start, end, '+' if strand == '+' else '-')
#
#             # Compare with corresponding reference sequence
#             if feature.featuretype == "miRNA_primary_transcript" and feature_id in pre_miRNA_ref:
#                 compare_sequences("miRNA_primary_transcript", extracted_seq, pre_miRNA_ref[feature_id], feature_id, output_handle, strand)
#             elif feature.featuretype == "miRNA" and feature_id in mature_star_ref:
#                 compare_sequences("miRNA", extracted_seq, mature_star_ref[feature_id], feature_id, output_handle, strand)
#
#     # Write summary
#     output_handle.write("Summary:\n")
#     output_handle.write(f"Total Matches: {match_count}\n")
#     output_handle.write(f"Total Mismatches: {mismatch_count}\n")
#
# print("Processing complete. Results saved in:", output_file)

import gffutils
from Bio import SeqIO
import os
import pandas as pd

# File paths
oscar_dir = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/miRge/"
genome_fasta = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/Genome/refs/Hmia_ref/Hmia.030120.fasta"
gff_file = oscar_dir + "miRNA_250203_reformatted.gff3"
pre_miRNA_fasta = oscar_dir + "pre_miR_1050_no_pre_in_seqid.fa"
mature_star_fasta = oscar_dir + "combined_mature_star_1050.fa"
output_file = oscar_dir + "genome_to_fasta_comparison_new.csv"
db_file = oscar_dir + "miRNA_250203_reformatted.db"

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

