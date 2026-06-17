"""Optional BLAST result loading and merge for intersections tables."""

import pandas as pd

BLAST_COLUMNS = [
    "query_accession", "subject_accession", "%_identical_matches", "alignment_length",
    "mismatches", "gap_openings", "query_start", "query_end",
    "subject_start", "subject_end", "e_value", "bitscore",
]


def load_blast_results(path, tool="mirdeep"):
    """Load BLAST tab output, dedupe by query, infer strand, extract index."""
    blast_orig = pd.read_csv(path, sep="\t", names=BLAST_COLUMNS)
    blast_orig = blast_orig.drop_duplicates(subset=["query_accession"])
    mask = blast_orig["subject_start"] < blast_orig["subject_end"]
    blast_orig.loc[mask, "strand"] = "+"
    blast_orig["strand"].fillna("-", inplace=True)
    blast = blast_orig.drop(
        ["%_identical_matches", "mismatches", "gap_openings", "subject_start", "subject_end", "bitscore"],
        axis=1,
    )
    index_pos = 4 if tool == "mirdeep" else 3
    blast["index"] = blast["query_accession"].str.split("|")
    blast["index"] = blast["index"].apply(lambda x: x[index_pos])
    return blast_orig, blast


def merge_blast(intersections_table, blast_df, description_col, index_pos=3):
    """Merge BLAST columns into an intersections table on extracted index."""
    table = intersections_table.copy()
    table["index"] = table[description_col].str.split(";")
    table["index"] = table["index"].apply(lambda x: x[index_pos])
    return pd.merge(table, blast_df, on="index", how="left")
