import argparse
import os

import gffutils
import pandas as pd
from Bio import SeqIO


def load_reference_sequences(fasta_file):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}


def extract_sequence(genome, seq_id, start, end, strand):
    sequence = genome[seq_id].seq[start - 1 : end]
    return str(sequence.reverse_complement()) if strand == "-" else str(sequence)


def compare_sequences(results, feature_type, extracted_seq, ref_seq, feature_id, strand):
    match_status = "Match" if extracted_seq == ref_seq else "Mismatch"
    results.append([feature_id, feature_type, extracted_seq, ref_seq, strand, match_status])


def verify_mirge_gff(db, genome, pre_miRNA_ref, mature_star_ref, results):
    for feature in db.all_features():
        if feature.featuretype not in ["miRNA_primary_transcript", "miRNA"]:
            continue
        feature_id = feature.attributes.get("ID", ["Unknown"])[0]
        extracted = extract_sequence(genome, feature.seqid, feature.start, feature.end, feature.strand)
        if feature.featuretype == "miRNA_primary_transcript" and feature_id in pre_miRNA_ref:
            compare_sequences(results, "miRNA_primary_transcript", extracted, pre_miRNA_ref[feature_id], feature_id, feature.strand)
        elif feature.featuretype == "miRNA" and feature_id in mature_star_ref:
            compare_sequences(results, "miRNA", extracted, mature_star_ref[feature_id], feature_id, feature.strand)


def verify_discovery_gff(db, genome, hairpin_by_name, mature_star_ref, results):
    for feature in db.all_features():
        feature_id = feature.attributes.get("ID", ["Unknown"])[0]
        extracted = extract_sequence(genome, feature.seqid, feature.start, feature.end, feature.strand)
        if feature.featuretype == "pre_miRNA":
            ref = hairpin_by_name.get(feature_id)
            if ref is not None:
                compare_sequences(results, "pre_miRNA", extracted, ref, feature_id, feature.strand)
        elif feature.featuretype == "miRNA":
            ref = mature_star_ref.get(feature_id)
            if ref is not None:
                compare_sequences(results, "miRNA", extracted, ref, feature_id, feature.strand)


def main():
    parser = argparse.ArgumentParser(description="Compare genome coordinates to reference FASTA sequences.")
    parser.add_argument("--species", required=True)
    parser.add_argument("--dir", required=True, help="Working directory for GFF/FASTA/db files")
    parser.add_argument("--genome_fasta", required=True)
    parser.add_argument("--gff", required=True, help="GFF filename within --dir")
    parser.add_argument("--premir", help="Hairpin/pre-miRNA FASTA filename within --dir")
    parser.add_argument("--mature_star", help="Mature/star combined FASTA filename within --dir")
    parser.add_argument("--hairpin", help="Discovery-mode hairpin FASTA (defaults to --premir)")
    parser.add_argument("--mature", help="Discovery-mode mature FASTA filename within --dir")
    parser.add_argument("--star", help="Discovery-mode star FASTA filename within --dir")
    parser.add_argument("--output", required=True)
    parser.add_argument("--db", help="gffutils DB filename within --dir (default: derived from GFF name)")
    parser.add_argument(
        "--hairpin-table",
        help="TSV/CSV with name and hairpinSeq columns for pre_miRNA verification (discovery mode)",
    )
    parser.add_argument(
        "--mode",
        choices=["mirge", "discovery"],
        default="mirge",
        help="mirge: miRNA_primary_transcript/miRNA; discovery: pre_miRNA/miRNA from tool GFF",
    )
    args = parser.parse_args()

    gff_file = os.path.join(args.dir, args.gff)
    output_file = os.path.join(args.dir, args.output)
    db_file = os.path.join(args.dir, args.db or (os.path.splitext(os.path.basename(args.gff))[0] + ".db"))

    genome = SeqIO.to_dict(SeqIO.parse(args.genome_fasta, "fasta"))
    results = []

    if not os.path.exists(db_file):
        print("Creating GFF database...")
        gffutils.create_db(
            gff_file,
            db_file,
            force=True,
            keep_order=True,
            merge_strategy="merge",
            sort_attribute_values=True,
        )
    db = gffutils.FeatureDB(db_file)

    if args.mode == "discovery":
        mature_file = os.path.join(args.dir, args.mature or args.mature_star)
        star_file = os.path.join(args.dir, args.star) if args.star else None
        mature_ref = load_reference_sequences(mature_file)
        star_ref = load_reference_sequences(star_file) if star_file and os.path.exists(star_file) else {}
        mature_star_ref = {**mature_ref, **star_ref}
        hairpin_by_name = {}
        if args.hairpin_table:
            hairpin_df = pd.read_csv(args.hairpin_table, sep="\t")
            if "name" in hairpin_df.columns and "hairpinSeq" in hairpin_df.columns:
                hairpin_by_name = dict(zip(hairpin_df["name"], hairpin_df["hairpinSeq"]))
        verify_discovery_gff(db, genome, hairpin_by_name, mature_star_ref, results)
    else:
        pre_miRNA_fasta = os.path.join(args.dir, args.premir)
        mature_star_fasta = os.path.join(args.dir, args.mature_star)
        pre_miRNA_ref = load_reference_sequences(pre_miRNA_fasta)
        mature_star_ref = load_reference_sequences(mature_star_fasta)
        verify_mirge_gff(db, genome, pre_miRNA_ref, mature_star_ref, results)

    df = pd.DataFrame(
        results,
        columns=["Feature ID", "Feature Type", "Extracted Sequence", "Reference Sequence", "Strand", "Match Status"],
    )
    df.to_csv(output_file, index=False)
    mismatches = (df["Match Status"] == "Mismatch").sum() if not df.empty else 0
    print(f"Processing complete. {len(df)} features checked, {mismatches} mismatches.")
    print("Results saved in:", output_file)


if __name__ == "__main__":
    main()
