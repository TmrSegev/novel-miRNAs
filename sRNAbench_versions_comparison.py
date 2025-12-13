#!/usr/bin/python

import sys
import pandas as pd
import argparse
import os
from Bio import SeqIO
import re

def load_fasta_sequences(fasta_file):
    """
    Load sequences from a FASTA file into a dictionary.
    Uses index for matching since IDs may differ between FASTA and GFF3.
    
    :param fasta_file: Path to FASTA file
    :return: Dictionary mapping index to sequence string
    """
    sequences = {}
    if not os.path.exists(fasta_file):
        print(f"Warning: FASTA file not found: {fasta_file}", file=sys.stderr)
        return sequences
    
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # FASTA header format: "new-mir-novel9-5p_1|m|1|index=1|MIR-71"
            # Extract index value for matching
            header = record.id  # BioPython record.id is the part before first space, or whole line if no space
            
            index = None
            # Look for index=X pattern
            index_match = re.search(r'index=(\d+)', header)
            if index_match:
                index = index_match.group(1)
            
            if index is not None:
                sequences[index] = str(record.seq)
            else:
                # Fallback: try to extract ID if no index found
                if 'ID=' in header:
                    id_match = re.search(r'ID=([^|;]+)', header)
                    if id_match:
                        sequences[id_match.group(1).strip()] = str(record.seq)
                else:
                    seq_id = header.split('|')[0].strip() if '|' in header else header.split(';')[0].strip()
                    if ' ' in seq_id:
                        seq_id = seq_id.split(' ')[0].strip()
                    if seq_id:
                        sequences[seq_id] = str(record.seq)
    except Exception as e:
        print(f"Error reading FASTA file {fasta_file}: {e}", file=sys.stderr)
    
    return sequences

def extract_index_from_attributes(attributes):
    """
    Extract index from GFF3 attributes column for matching with FASTA sequences.
    
    :param attributes: Attributes string (e.g., "ID=new-mir-novel9_1;RC_m=11;RC_s=6;index=1;MIR-71;novel;2")
    :return: Index string or None (e.g., "1")
    """
    if pd.isna(attributes):
        return None
    match = re.search(r'index=(\d+)', str(attributes))
    if match:
        return match.group(1)
    return None

def extract_id_from_attributes(attributes):
    """
    Extract ID from GFF3 attributes column.
    
    :param attributes: Attributes string (e.g., "ID=new-mir-novel9_1;RC_m=11;RC_s=6;index=1;MIR-71;novel;2")
    :return: ID string or None (e.g., "new-mir-novel9_1")
    """
    if pd.isna(attributes):
        return None
    match = re.search(r'ID=([^;]+)', str(attributes))
    if match:
        return match.group(1).strip()
    return None

def count_intersections(intersections_file, debug=False, output_non_intersected=None, fasta_file=None):
    """
    Count how many miRNA candidates have an intersection in the BED file.
    
    :param intersections_file: Path to the bedtools intersection BED file
    :param debug: If True, print debug information
    :param output_non_intersected: Optional path to save non-intersected candidates
    :param fasta_file: Optional path to FASTA file to add sequences to output
    :return: Dictionary with counts and statistics
    """
    # Read the BED file
    # BED format from bedtools intersect has:
    # Columns 1-9: First file (v2)
    # Columns 10-18: Second file (v1) - intersection data
    column_names = [
        'Chr_v2', '.1', 'type_v2', 'Start_v2', 'End_v2', '.2', 'Strand_v2', '.3', 'Description_v2',
        'Chr_v1', '.4', 'type_v1', 'Start_v1', 'End_v1', '.5', 'Strand_v1', '.6', 'Description_v1'
    ]
    
    try:
        # Read without column names first to see actual structure
        df_raw = pd.read_csv(intersections_file, sep='\t', header=None)
    except Exception as e:
        print(f"Error reading file {intersections_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Check if file is empty
    if df_raw.empty:
        print("Warning: Intersections file is empty.")
        return {
            'total_candidates': 0,
            'with_intersection': 0,
            'without_intersection': 0,
            'intersection_rate': 0.0
        }
    
    # Determine if there's an intersection
    # When there's NO intersection, bedtools outputs: columns 1-9 from first file, then . . . -1 -1 . . . . 0
    # When there IS an intersection, bedtools outputs: columns 1-9 from first file, then columns 1-9 from second file
    # So we check column 10 (index 9) - if it's a valid scaffold/chromosome name (not '.', '-1', '0', or empty), there's an intersection
    
    # Column 10 is at index 9 (0-based)
    if df_raw.shape[1] <= 9:
        # File doesn't have enough columns, no intersections possible
        col_10 = pd.Series(['.'] * len(df_raw))
    else:
        col_10 = df_raw.iloc[:, 9].astype(str).str.strip()
    
    if debug:
        print(f"\nDebug: Total columns in file: {df_raw.shape[1]}", file=sys.stderr)
        print(f"Debug: First 5 values in column 10:", file=sys.stderr)
        print(col_10.head().tolist(), file=sys.stderr)
        print(f"Debug: Unique values in column 10 (first 20):", file=sys.stderr)
        print(col_10.unique()[:20], file=sys.stderr)
    
    # An intersection exists if column 10 contains a scaffold/chromosome name
    # Valid intersections have scaffold names (like "scaffold_94"), not '.', '-1', '0', or empty
    # Simple check: if it's not one of the invalid markers, it's an intersection
    invalid_markers = {'.', '-1', '0', '', 'nan', 'NaN', 'None'}
    has_intersection = ~col_10.isin(invalid_markers) & col_10.notna()
    
    total_candidates = len(df_raw)
    with_intersection = has_intersection.sum()
    without_intersection = total_candidates - with_intersection
    intersection_rate = (with_intersection / total_candidates * 100) if total_candidates > 0 else 0.0
    
    # Save non-intersected candidates if requested
    non_intersected_saved = False
    if output_non_intersected and without_intersection > 0:
        # Extract non-intersected candidates (columns 1-9 from first file)
        non_intersected = df_raw[~has_intersection].iloc[:, :9].copy()
        
        # Add column names for the first 9 columns (GFF3 format)
        non_intersected.columns = [
            'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'
        ]
        
        # Load sequences from FASTA file if provided
        sequences_dict = {}
        if fasta_file and os.path.exists(fasta_file):
            sequences_dict = load_fasta_sequences(fasta_file)
            if debug:
                print(f"Debug: Loaded {len(sequences_dict)} sequences from FASTA file: {fasta_file}", file=sys.stderr)
                if len(sequences_dict) > 0:
                    # Show sample indices from FASTA (keys are now indices, not IDs)
                    sample_indices = list(sequences_dict.keys())[:10]
                    print(f"Debug: Sample FASTA indices (keys): {sample_indices}", file=sys.stderr)
        elif fasta_file:
            print(f"Warning: FASTA file not found: {fasta_file}", file=sys.stderr)
        
        # Add sequence column
        if sequences_dict:
            # Extract index from attributes and match with sequences (index is more reliable than ID)
            non_intersected['sequence'] = non_intersected['attributes'].apply(
                lambda attr: sequences_dict.get(extract_index_from_attributes(attr), '')
            )
            matched_count = (non_intersected['sequence'] != '').sum()
            if debug:
                # Show sample indices from attributes
                sample_attrs = non_intersected['attributes'].head(5).tolist()
                sample_extracted_indices = [extract_index_from_attributes(attr) for attr in sample_attrs]
                sample_extracted_ids = [extract_id_from_attributes(attr) for attr in sample_attrs]
                print(f"Debug: Sample extracted indices from attributes: {sample_extracted_indices}", file=sys.stderr)
                print(f"Debug: Sample extracted IDs from attributes: {sample_extracted_ids}", file=sys.stderr)
                print(f"Debug: Sample indices in FASTA dict: {list(sequences_dict.keys())[:10]}", file=sys.stderr)
                print(f"Debug: Matched {matched_count} out of {len(non_intersected)} candidates with sequences", file=sys.stderr)
            if matched_count == 0:
                print(f"Warning: No sequences matched! Check if index format in FASTA headers matches GFF3 attributes.", file=sys.stderr)
        else:
            # Add empty sequence column if FASTA not provided or not found
            non_intersected['sequence'] = ''
            if fasta_file:
                print(f"Warning: No sequences loaded. FASTA file may be empty or not found.", file=sys.stderr)
        
        # Save to CSV
        non_intersected.to_csv(output_non_intersected, sep=',', index=False)
        non_intersected_saved = True
        if debug:
            print(f"Debug: Saved {without_intersection} non-intersected candidates to {output_non_intersected}", file=sys.stderr)
    
    return {
        'total_candidates': total_candidates,
        'with_intersection': with_intersection,
        'without_intersection': without_intersection,
        'intersection_rate': intersection_rate,
        'non_intersected_file': output_non_intersected if non_intersected_saved else None
    }


def main():
    parser = argparse.ArgumentParser(
        description='Count miRNA candidates with intersections in sRNAbench versions comparison BED file'
    )
    parser.add_argument(
        '--intersections-file',
        required=True,
        help='Path to the bedtools intersection BED file (e.g., sRNAbench_intersect_v2_v1.bed)'
    )
    parser.add_argument(
        '--output',
        help='Optional: Path to save detailed results as CSV (default: print to stdout)'
    )
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Print debug information about column values'
    )
    parser.add_argument(
        '--output-non-intersected',
        help='Path to save candidates without intersections (CSV format, comma-separated)'
    )
    parser.add_argument(
        '--fasta-file',
        help='Path to FASTA file to add sequences to output (e.g., Hofstenia_sRNAbench_pre_only_v2.fasta). If not provided, will auto-detect from intersection filename.'
    )
    parser.add_argument(
        '--fasta-version',
        choices=['v1', 'v2'],
        help='Manually specify FASTA file version (v1 or v2). If not provided, will auto-detect from intersection filename (e.g., sRNAbench_intersect_v2_v1.bed -> v2).'
    )
    
    args = parser.parse_args()
    
    # Auto-detect FASTA file version from intersection file name
    # Format: sRNAbench_intersect_v2_v1.bed means v2 is query (first file), so use v2 FASTA
    # Format: sRNAbench_intersect_v1_v2.bed means v1 is query (first file), so use v1 FASTA
    fasta_file = args.fasta_file
    auto_detected_version = None
    
    if not fasta_file:
        # Try to auto-detect version from intersection file name
        intersection_basename = os.path.basename(args.intersections_file)
        # Look for pattern: *_vX_vY.bed or *_vX_vY.*
        version_match = re.search(r'_v(\d+)_v(\d+)', intersection_basename)
        if version_match:
            # First version in filename is the query (first file in intersection)
            auto_detected_version = f"v{version_match.group(1)}"
            if args.debug:
                print(f"Debug: Auto-detected FASTA version: {auto_detected_version} from filename", file=sys.stderr)
        elif args.fasta_version:
            # Use explicitly provided version
            auto_detected_version = args.fasta_version
    
    # Build FASTA file path if version is known
    if auto_detected_version and not fasta_file:
        # Try to find FASTA file in same directory as intersection file
        intersections_dir = os.path.dirname(os.path.abspath(args.intersections_file))
        fasta_file = os.path.join(intersections_dir, f'Hofstenia_sRNAbench_pre_only_{auto_detected_version}.fasta')
        if not os.path.exists(fasta_file):
            # Try alternative path
            scripts_dir = '/groups/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia_newGenome/scripts/'
            fasta_file = os.path.join(scripts_dir, f'Hofstenia_sRNAbench_pre_only_{auto_detected_version}.fasta')
        if not os.path.exists(fasta_file) and args.debug:
            print(f"Debug: Could not find FASTA file: {fasta_file}", file=sys.stderr)
    
    # Count intersections
    results = count_intersections(args.intersections_file, debug=args.debug, 
                                   output_non_intersected=args.output_non_intersected,
                                   fasta_file=fasta_file)
    
    # Print results
    print("=" * 60)
    print("sRNAbench Versions Comparison Results")
    print("=" * 60)
    print(f"Total miRNA candidates: {results['total_candidates']}")
    print(f"Candidates with intersection: {results['with_intersection']}")
    print(f"Candidates without intersection: {results['without_intersection']}")
    print(f"Intersection rate: {results['intersection_rate']:.2f}%")
    print("=" * 60)
    
    # Report if non-intersected candidates were saved
    if results['non_intersected_file']:
        print(f"\nNon-intersected candidates saved to: {results['non_intersected_file']}")
        print(f"  ({results['without_intersection']} candidates)")
        if fasta_file:
            if os.path.exists(fasta_file):
                print(f"  FASTA file used: {fasta_file}")
                # Check if sequences were actually added by reading the output file
                try:
                    output_df = pd.read_csv(results['non_intersected_file'], sep=',')
                    if 'sequence' in output_df.columns:
                        seq_count = (output_df['sequence'] != '').sum()
                        print(f"  Sequences matched: {seq_count} out of {len(output_df)} candidates")
                except:
                    pass
            else:
                print(f"  Warning: FASTA file not found: {fasta_file}")
    
    # Save detailed results if output path is provided
    if args.output:
        results_df = pd.DataFrame([results])
        results_df.to_csv(args.output, index=False)
        print(f"\nDetailed results saved to: {args.output}")


if __name__ == '__main__':
    main()

