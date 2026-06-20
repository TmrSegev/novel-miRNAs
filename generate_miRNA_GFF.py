import pandas as pd
from pathlib import Path
import os
import argparse

from pipeline_config import get_species_config

def locate_subseq(full_seq, sub_seq, strand, start_offset, end_offset): ## pass end_offset
    idx = full_seq.find(sub_seq)
    if idx == -1:
        raise ValueError(f"Could not find subsequence: {sub_seq} in {full_seq}")
    
    if strand == '+':
        gff_start = start_offset + idx
        gff_end = gff_start + len(sub_seq) - 1        
    else: ## new
        gff_end = end_offset - idx
        gff_start = gff_end - len(sub_seq) + 1
        
    return gff_start + 1, gff_end + 1  # 1-based GFF3 format

def create_gff_rows(row):
    rows = []
    seqname = row['Chr']
    source = "."
    strand = row['Strand']
    start = int(row['Start'])
    end = int(row['End'])
    hairpin = str(row.get('hairpinSeq', '')).strip()
    five_p_seq = str(row.get('5pseq', '')).strip()
    three_p_seq = str(row.get('3pseq', '')).strip()
    mirna_name = row['Description'].replace(';', '|')
    mature_type = row['mature']  # 5p or 3p

    if not all([hairpin, five_p_seq, three_p_seq]):
        raise ValueError("One or more sequences are missing")

    # Assign mature and star sequences
    if mature_type == "5p":
        mature_seq = five_p_seq
        star_seq = three_p_seq
    else:
        mature_seq = three_p_seq
        star_seq = five_p_seq

    # Primary transcript
    attr = f"ID={mirna_name};Alias={mirna_name};Name={mirna_name}"
    rows.append([seqname, source, "miRNA_primary_transcript", start, end, ".", strand, ".", attr])

    # Mature miRNA
    m_start, m_end = locate_subseq(hairpin, mature_seq, strand, start - 1, end - 1) ## pass the end
    m_id = f"{mirna_name}-{mature_type}"
    m_attr = f"ID={m_id};Alias={m_id};Name={m_id};Derives_from={mirna_name}"
    #rows.append([seqname, source, "miRNA", m_start, m_end, ".", strand, ".", m_attr])

    # Star miRNA
    s_type = "3p" if mature_type == "5p" else "5p"
    s_start, s_end = locate_subseq(hairpin, star_seq, strand, start - 1, end - 1) ## pass the end
    s_id = f"{mirna_name}-{s_type}"
    s_attr = f"ID={s_id};Alias={s_id};Name={s_id};Derives_from={mirna_name}"

    ### print first the one that appears first
    if m_start < s_start:
       rows.append([seqname, source, "miRNA", m_start, m_end, ".", strand, ".", m_attr])
       rows.append([seqname, source, "miRNA", s_start, s_end, ".", strand, ".", s_attr])   
    else:
       rows.append([seqname, source, "miRNA", s_start, s_end, ".", strand, ".", s_attr])   
       rows.append([seqname, source, "miRNA", m_start, m_end, ".", strand, ".", m_attr])
        
    return rows


# --- Argument Parsing ---
# Set up a parser to read command-line arguments
parser = argparse.ArgumentParser(
    description="Combines 'all_candidates_mature.fasta' and 'all_candidates_star.fasta' from a base directory into a single FASTA file. It adds ';mature' and ';star' suffixes to the respective sequence headers."
)
# Add a required positional argument for the base directory path
parser.add_argument(
    "--species",
    help="Species name (or Species_newGenome alias for alternate assembly).",
)
parser.add_argument("--variant", help='Genome variant, e.g. "new_genome"')
parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
# Parse the arguments provided by the user
args = parser.parse_args()
SPECIES = args.species
try:
    cfg = get_species_config(SPECIES, args.base_path, variant=args.variant)
except ValueError:
    cfg = None

output_track = cfg["variant_track"] if cfg and cfg.get("variant") == "new_genome" else SPECIES
input_excel = Path(
    f"/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/"
    f"all_remaining_after_ziv_{output_track}.xlsx"
)
sheet_name = cfg["mirge_input_sheet"] if cfg else "(D) Structural Features"

if cfg and cfg.get("variant") == "new_genome":
    mirge_root = Path(cfg["variant_root_dir"])
elif cfg:
    mirge_root = Path(cfg["base_path"], cfg["species"])
else:
    mirge_root = Path(f"/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/{SPECIES}")

if cfg and cfg.get("ziv_profile") == "unfiltered_only":
    output_gff = mirge_root / "miRge_after_Ziv" / "miRNA_candidates.gff3"
else:
    output_gff = mirge_root / "miRge" / "miRNA_candidates.gff3"

# Load Excel data
df = pd.read_excel(input_excel, sheet_name=sheet_name, dtype=str)

# Build GFF rows
gff_rows = []
for _, row in df.iterrows():
    try:
        gff_rows.extend(create_gff_rows(row))
    except ValueError as e:
        print(f"Warning: {e} (miRNA: {row['Description']})")

# Save as GFF3
gff_df = pd.DataFrame(gff_rows, columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])
with open(output_gff, "w") as f:
    gff_df.to_csv(f, sep="\t", header=False, index=False, line_terminator="\n")

print(f"Saved GFF3 to: {output_gff}")
