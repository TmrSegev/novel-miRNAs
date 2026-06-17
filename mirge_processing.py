import os
import pandas as pd
import re
import argparse
from pathlib import Path

from pipeline_config import SPECIES_CONFIG, get_species_config

# --- Argument Parsing ---
# Set up a parser to read command-line arguments
parser = argparse.ArgumentParser(
    description=""
)

# Add a required positional argument for the base directory path
parser.add_argument("--dir", help="")
parser.add_argument("--species", required=True, help="Species name")
parser.add_argument("--variant", help='Genome variant, e.g. "new_genome"')
parser.add_argument("--base-path", dest="base_path", help="Charles_seq base path")
# Add optional flag to toggle *_m18 naming
parser.add_argument("--m18", action="store_true", help="If set, use *_m18 paths and filenames")

# Parse the arguments provided by the user
args = parser.parse_args()
m18_suffix = "_m18" if args.m18 else ""
dir = args.dir
SPECIES = args.species
cfg = None
if SPECIES in SPECIES_CONFIG:
    cfg = get_species_config(SPECIES, args.base_path, variant=args.variant)
elif SPECIES == "Hofstenia_newGenome":
    cfg = get_species_config("Hofstenia", args.base_path, variant="new_genome")

sheet_name = cfg["mirge_input_sheet"] if cfg else "(D) Structural Features"

input_excel = Path(
    f"/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_{SPECIES}.xlsx")

print(f"Loading maturity info from: {input_excel} | Sheet: {sheet_name}")
# Load the maturity data from the specified sheet in the Excel file
maturity_df = pd.read_excel(input_excel, sheet_name=sheet_name)

# Keep mature and (if present) Strand
cols_to_keep = ['Description', 'mature'] + (['Strand'] if 'Strand' in maturity_df.columns else [])
maturity_df = maturity_df[cols_to_keep].copy()

# Build maps (by precursor / Description)
maturity_map = pd.Series(
    maturity_df.dropna(subset=['mature'])['mature'].values,
    index=maturity_df.dropna(subset=['mature'])['Description'].str.strip()
).to_dict()

strand_map = {}
if 'Strand' in maturity_df.columns:
    strand_map = pd.Series(
        maturity_df.dropna(subset=['Strand'])['Strand'].values,
        index=maturity_df.dropna(subset=['Strand'])['Description'].str.strip()
    ).to_dict()


# Define the base directory containing all libraries
# base_dir = '/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/results/'
base_dir = os.path.join(dir, f'results{m18_suffix}/')

# Variables to track the number of miRNAs missing 3p or 5p across all libraries
no_3p_or_5p_mirs_overall = 0

# List to store all paired data from different libraries
all_paired_data = []
all_data = []

# Iterate through each library folder
for library_folder in os.listdir(base_dir):
    library_path = os.path.join(base_dir, library_folder)
    gff_file_path = os.path.join(library_path, 'sample_miRge3.gff')

    # Skip if "sample_miRge3.gff" is not found in the directory
    if not os.path.isfile(gff_file_path):
        continue

    no_3p_or_5p_mirs_library = 0

    # Extract library name from the third line
    library_name = None
    with open(gff_file_path, 'r') as f:
        for line in f:
            if line.startswith("## COLDATA:"):
                match = re.search(r'## COLDATA:\s*([\w-]+)', line)
                if match:
                    library_name = match.group(1).split('.')[0]  # Extract before the dot
                    break

    if library_name is None:
        raise ValueError(f"Library name not found in file: {gff_file_path}")

    data = []

    with open(gff_file_path, 'r') as f:
        for _ in range(4):  # Skip the first 4 lines
            next(f)
        for line in f:
            parts = line.strip().split('\t')
            attributes = dict(item.split('=', 1) for item in parts[-1].split(';'))
            attributes = {key.lstrip(): value for key, value in attributes.items()}
            data.append([
                parts[0],  # miRNA name
                parts[2],  # Type (ref_miRNA or isomiR)
                int(attributes['Hits']),
                attributes.get('Variant', 'NA'),
                attributes.get('Read', 'NA'),  # Extract Sequence
                len(attributes.get('Read', '')) if 'Read' in attributes else 'NA'  # Length of sequence
            ])

    df = pd.DataFrame(data, columns=['miRNA', 'Type', 'Hits', 'Variant', 'Sequence', 'Length'])

    # Add library name as a new column
    df['Library'] = library_name


    def process_mirna_group(group):
        total_hits = group['Hits'].sum()
        total_hits_filtered = group[~group['Variant'].str.contains('snv|iso_add5p|iso_add3p', na=False)]['Hits'].sum()
        ref_hits = group[group['Type'] == 'ref_miRNA']['Hits'].sum()

        higher_isomir = 'NA'
        higher_isomir_seq = 'NA'
        higher_isomir_length = 'NA'
        higher_isomir_count = 'NA'

        # Filter out isomiRs containing snv, iso_add5p, or iso_add3p before finding the highest
        isomir_group = group[(group['Type'] == 'isomiR') & 
                             ~group['Variant'].str.contains('snv|iso_add5p|iso_add3p', na=False)]
        if not isomir_group.empty:
            max_isomir_row = isomir_group.loc[isomir_group['Hits'].idxmax()]
            if ref_hits == 0 or max_isomir_row['Hits'] > ref_hits:
                higher_isomir = max_isomir_row['Variant']
                higher_isomir_seq = max_isomir_row['Sequence']
                higher_isomir_length = max_isomir_row['Length']
                higher_isomir_count = max_isomir_row['Hits']

        iso_5p_hits = group[group['Variant'].str.contains(r'iso_5p:', na=False) &
                            ~group['Variant'].str.contains('snv|iso_add5p|iso_add3p', na=False)]['Hits'].sum()
        iso_3p_hits = group[group['Variant'].str.contains(r'iso_3p', na=False) &
                            ~group['Variant'].str.contains('snv|iso_add5p|iso_add3p', na=False)]['Hits'].sum()
        iso_add5p_hits = group[group['Variant'].str.contains('iso_add5p', na=False)]['Hits'].sum()
        iso_add3p_hits = group[group['Variant'].str.contains('iso_add3p', na=False)]['Hits'].sum()

        percent_5p = (iso_5p_hits / total_hits_filtered) * 100 if total_hits_filtered > 0 else 0
        percent_3p = (iso_3p_hits / total_hits_filtered) * 100 if total_hits_filtered > 0 else 0

        return [
            group['miRNA'].iloc[0], total_hits, total_hits_filtered, ref_hits, higher_isomir,
            higher_isomir_seq, higher_isomir_length, higher_isomir_count,
            iso_5p_hits, iso_3p_hits, iso_add5p_hits, iso_add3p_hits, percent_5p, percent_3p, group['Library'].iloc[0]
        ]


    report_data = [process_mirna_group(group) for _, group in df.groupby('miRNA')]

    columns = [
        'miRNA', 'Total Hits', 'Total Hits Filtered', 'ref_miRNA Hits', 'Higher isomiR Variant',
        'Higher isomiR Sequence', 'Higher isomiR Length', 'Higher isomiR Count',
        '5p isomiRs Count', '3p isomiRs Count', '5add isomiRs Count', '3add isomiRs Count',
        '% of 5p isomiRs / Total Filtered', '% of 3p isomiRs / Total Filtered', 'Library'
    ]

    report_df = pd.DataFrame(report_data, columns=columns)
    all_data.append(report_df)

# Concatenate all paired data into a final dataframe
final_df = pd.concat(all_data, ignore_index=True)
final_df = final_df.drop_duplicates(subset=['miRNA', 'Library'], keep='last')

# Sort by pre-miRNA and Library
final_df.sort_values(by=['miRNA', 'Library'], inplace=True)


# --- NEW: Add 'mature' Column ---
# Define a function to determine if a miRNA is 'mature' or 'Star'
def get_maturity_status(row, mapping):
    full_mirna_name = row['miRNA']
    # Split the name into precursor and arm (e.g., '...-3p' -> ['...', '3p'])
    parts = full_mirna_name.rsplit('-', 1)

    # Check if the split was successful and gave a valid arm
    if len(parts) < 2 or parts[1] not in ['5p', '3p']:
        return 'Unknown'  # Cannot determine arm

    precursor_name, current_arm = parts[0], parts[1]

    # Look up the designated mature arm from the map created earlier
    designated_mature_arm = mapping.get(precursor_name)

    if designated_mature_arm is None:
        return 'Unknown'  # Precursor not found in Excel file

    # Compare the miRNA's arm with the designated mature arm
    return 'mature' if current_arm == designated_mature_arm else 'star'


def get_strand(row, mapping):
    full_mirna_name = row['miRNA']
    parts = full_mirna_name.rsplit('-', 1)
    if len(parts) < 2 or parts[1] not in ['5p', '3p']:
        return 'Unknown'
    precursor_name = parts[0]
    return mapping.get(precursor_name, 'Unknown')



# Apply the function to each row to create the new 'mature' column
final_df['mature'] = final_df.apply(lambda row: get_maturity_status(row, maturity_map), axis=1)
# Add Strand column (if Strand info exists in the Excel)
final_df['Strand'] = final_df.apply(lambda row: get_strand(row, strand_map), axis=1) if strand_map else 'Unknown'


cols = final_df.columns.tolist()
# Put 'mature' right after 'miRNA'
if 'mature' in cols:
    cols.insert(1, cols.pop(cols.index('mature')))
# Put 'Strand' right after 'mature'
if 'Strand' in cols:
    cols.insert(2, cols.pop(cols.index('Strand')))
final_df = final_df[cols]

# --- END NEW SECTION ---

# Save the final concatenated and sorted paired dataframe
output_dir = os.path.join(dir, f"processing_results{m18_suffix}/")
os.makedirs(output_dir, exist_ok=True)
final_output_filename = os.path.join(output_dir, f'final_miRNA_isomiR_report{m18_suffix}.csv')
final_df.to_csv(final_output_filename, index=False)

print("FINISH processing final_df")

# --- NEW: Enhance input_excel with isomiR percentages ---
print("Starting enhancement of input Excel file...")

# Load the original Excel file
original_df = pd.read_excel(input_excel, sheet_name=sheet_name)
print(f"Original Excel has {len(original_df)} rows")

# Get unique libraries from final_df
libraries = sorted(final_df['Library'].unique())
print(f"Found {len(libraries)} libraries: {libraries}")


# Function to extract precursor name from miRNA name
def extract_precursor_name(mirna_name):
    """Extract precursor name from miRNA name (remove -5p or -3p suffix)"""
    return mirna_name.rsplit('-', 1)[0] if mirna_name.endswith(('-5p', '-3p')) else mirna_name


# Create pivot tables for 5p and 3p percentages
print("Creating pivot tables for percentages...")

# Pivot for 5p percentages
pivot_5p = final_df.pivot_table(
    index='miRNA',
    columns='Library',
    values='% of 5p isomiRs / Total Filtered',
    fill_value=0
)

# Pivot for 3p percentages
pivot_3p = final_df.pivot_table(
    index='miRNA',
    columns='Library',
    values='% of 3p isomiRs / Total Filtered',
    fill_value=0
)


# Function to match candidates and add percentage data
def add_percentage_columns(original_df, pivot_5p, pivot_3p, libraries):
    """Add percentage columns for each library to the original dataframe"""
    enhanced_df = original_df.copy()

    # Initialize new columns with NaN
    for lib in libraries:
        enhanced_df[f'{lib}_5p_percent'] = float('nan')
        enhanced_df[f'{lib}_3p_percent'] = float('nan')

    # For each row in original dataframe
    for idx, row in original_df.iterrows():
        candidate_name = row['Description'].strip()

        # Look for matching miRNAs in the pivot tables
        # Try both 5p and 3p versions
        mirna_5p = f"{candidate_name}-5p"
        mirna_3p = f"{candidate_name}-3p"

        # Add 5p percentages if available
        if mirna_5p in pivot_5p.index:
            for lib in libraries:
                if lib in pivot_5p.columns:
                    enhanced_df.loc[idx, f'{lib}_5p_percent'] = pivot_5p.loc[mirna_5p, lib]

        # Add 3p percentages if available
        if mirna_3p in pivot_3p.index:
            for lib in libraries:
                if lib in pivot_3p.columns:
                    enhanced_df.loc[idx, f'{lib}_3p_percent'] = pivot_3p.loc[mirna_3p, lib]

    return enhanced_df


# Enhance the original dataframe
enhanced_df = add_percentage_columns(original_df, pivot_5p, pivot_3p, libraries)

# --- REMOVE LEGACY _percent COLUMNS, THEY'RE REPLACED BY MATURE/STAR COLUMNS ---
legacy_cols = []
for lib in libraries:
    legacy_cols.extend([f'{lib}_5p_percent', f'{lib}_3p_percent'])

# Drop only columns that actually exist
cols_to_drop = [c for c in legacy_cols if c in enhanced_df.columns]
if cols_to_drop:
    enhanced_df.drop(columns=cols_to_drop, inplace=True)
    print(f"Removed legacy columns: {len(cols_to_drop)} columns")
else:
    print("No legacy *_percent columns found to remove.")

# Create pivot tables for Mature and Star
pivot_mature_5p = final_df[final_df['mature'] == 'mature'].pivot_table(
    index='miRNA', columns='Library', values='% of 5p isomiRs / Total Filtered', fill_value=0
)
pivot_mature_3p = final_df[final_df['mature'] == 'mature'].pivot_table(
    index='miRNA', columns='Library', values='% of 3p isomiRs / Total Filtered', fill_value=0
)
pivot_mature_total = final_df[final_df['mature'] == 'mature'].pivot_table(
    index='miRNA', columns='Library', values='Total Hits Filtered', fill_value=0
)

pivot_star_5p = final_df[final_df['mature'] == 'star'].pivot_table(
    index='miRNA', columns='Library', values='% of 5p isomiRs / Total Filtered', fill_value=0
)
pivot_star_3p = final_df[final_df['mature'] == 'star'].pivot_table(
    index='miRNA', columns='Library', values='% of 3p isomiRs / Total Filtered', fill_value=0
)
pivot_star_total = final_df[final_df['mature'] == 'star'].pivot_table(
    index='miRNA', columns='Library', values='Total Hits Filtered', fill_value=0
)

# Initialize dictionary to collect all new columns at once
new_columns_data = {}
for lib in libraries:
    new_columns_data[f'Mature-T-{lib}'] = []
    new_columns_data[f'Mature-%5p-{lib}'] = []
    new_columns_data[f'Mature-%3p-{lib}'] = []
    new_columns_data[f'Star-T-{lib}'] = []
    new_columns_data[f'Star-%5p-{lib}'] = []
    new_columns_data[f'Star-%3p-{lib}'] = []

# Populate the columns efficiently
for idx, row in enhanced_df.iterrows():
    candidate_name = row['Description'].strip()
    mirna_5p = f"{candidate_name}-5p"
    mirna_3p = f"{candidate_name}-3p"

    for lib in libraries:
        # Initialize with NaN
        mature_t = mature_5p = mature_3p = float('nan')
        star_t = star_5p = star_3p = float('nan')
        
        # Mature data
        if mirna_5p in pivot_mature_total.index and lib in pivot_mature_total.columns:
            mature_t = pivot_mature_total.loc[mirna_5p, lib]
            mature_5p = pivot_mature_5p.loc[mirna_5p, lib]
            mature_3p = pivot_mature_3p.loc[mirna_5p, lib]
        elif mirna_3p in pivot_mature_total.index and lib in pivot_mature_total.columns:
            mature_t = pivot_mature_total.loc[mirna_3p, lib]
            mature_5p = pivot_mature_5p.loc[mirna_3p, lib]
            mature_3p = pivot_mature_3p.loc[mirna_3p, lib]

        # Star data
        if mirna_5p in pivot_star_total.index and lib in pivot_star_total.columns:
            star_t = pivot_star_total.loc[mirna_5p, lib]
            star_5p = pivot_star_5p.loc[mirna_5p, lib]
            star_3p = pivot_star_3p.loc[mirna_5p, lib]
        elif mirna_3p in pivot_star_total.index and lib in pivot_star_total.columns:
            star_t = pivot_star_total.loc[mirna_3p, lib]
            star_5p = pivot_star_5p.loc[mirna_3p, lib]
            star_3p = pivot_star_3p.loc[mirna_3p, lib]
        
        # Append to lists
        new_columns_data[f'Mature-T-{lib}'].append(mature_t)
        new_columns_data[f'Mature-%5p-{lib}'].append(mature_5p)
        new_columns_data[f'Mature-%3p-{lib}'].append(mature_3p)
        new_columns_data[f'Star-T-{lib}'].append(star_t)
        new_columns_data[f'Star-%5p-{lib}'].append(star_5p)
        new_columns_data[f'Star-%3p-{lib}'].append(star_3p)

# Add all columns at once using concat (much faster than repeated column assignment)
new_columns_df = pd.DataFrame(new_columns_data, index=enhanced_df.index)
enhanced_df = pd.concat([enhanced_df, new_columns_df], axis=1)

# Initialize lists for Max columns
mature_t_max_list = []
mature_5p_max_list = []
mature_t_max_lib_list = []
star_t_max_list = []
star_5p_max_list = []
star_t_max_lib_list = []

# Find the library with max Mature-T and Star-T for each row
for idx, row in enhanced_df.iterrows():
    # Find max Mature-T across all libraries
    mature_t_values = {}
    mature_5p_values = {}
    for lib in libraries:
        mature_t_col = f'Mature-T-{lib}'
        mature_5p_col = f'Mature-%5p-{lib}'
        if mature_t_col in enhanced_df.columns:
            mature_t_val = row[mature_t_col]
            mature_5p_val = row[mature_5p_col] if mature_5p_col in enhanced_df.columns else float('nan')
            if pd.notna(mature_t_val) and mature_t_val > 0:
                mature_t_values[lib] = mature_t_val
                mature_5p_values[lib] = mature_5p_val
    
    if mature_t_values:
        max_mature_lib = max(mature_t_values, key=mature_t_values.get)
        mature_t_max_list.append(mature_t_values[max_mature_lib])
        mature_5p_max_list.append(mature_5p_values[max_mature_lib])
        mature_t_max_lib_list.append(max_mature_lib)
    else:
        mature_t_max_list.append(float('nan'))
        mature_5p_max_list.append(float('nan'))
        mature_t_max_lib_list.append('')
    
    # Find max Star-T across all libraries
    star_t_values = {}
    star_5p_values = {}
    for lib in libraries:
        star_t_col = f'Star-T-{lib}'
        star_5p_col = f'Star-%5p-{lib}'
        if star_t_col in enhanced_df.columns:
            star_t_val = row[star_t_col]
            star_5p_val = row[star_5p_col] if star_5p_col in enhanced_df.columns else float('nan')
            if pd.notna(star_t_val) and star_t_val > 0:
                star_t_values[lib] = star_t_val
                star_5p_values[lib] = star_5p_val
    
    if star_t_values:
        max_star_lib = max(star_t_values, key=star_t_values.get)
        star_t_max_list.append(star_t_values[max_star_lib])
        star_5p_max_list.append(star_5p_values[max_star_lib])
        star_t_max_lib_list.append(max_star_lib)
    else:
        star_t_max_list.append(float('nan'))
        star_5p_max_list.append(float('nan'))
        star_t_max_lib_list.append('')

# Add Max columns all at once
enhanced_df['Mature-T-Max'] = mature_t_max_list
enhanced_df['Mature-%5p-Max'] = mature_5p_max_list
enhanced_df['Mature-T-Max-Library'] = mature_t_max_lib_list
enhanced_df['Star-T-Max'] = star_t_max_list
enhanced_df['Star-%5p-Max'] = star_5p_max_list
enhanced_df['Star-T-Max-Library'] = star_t_max_lib_list


# --- REORDER: Group by column type instead of by library ---
new_cols_order = []

# Define the 6 column types
column_types = ['Mature-T', 'Mature-%5p', 'Mature-%3p', 'Star-T', 'Star-%5p', 'Star-%3p']

# For each column type, add all libraries
for col_type in column_types:
    for lib in libraries:
        new_cols_order.append(f'{col_type}-{lib}')

# Add Max columns after all library-specific columns
new_cols_order.extend(['Mature-T-Max', 'Mature-%5p-Max', 'Mature-T-Max-Library', 'Star-T-Max', 'Star-%5p-Max', 'Star-T-Max-Library'])

existing_cols = [c for c in enhanced_df.columns if c not in new_cols_order]
enhanced_df = enhanced_df[existing_cols + new_cols_order]


# Save the enhanced dataframe
candidates_dir = f"/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral/miRNAs/{SPECIES}/"
output_excel_path = os.path.join(candidates_dir, f'final_candidates_v1{m18_suffix}.xlsx')
enhanced_df.to_excel(output_excel_path, index=False)

print(f"Enhanced Excel saved to: {output_excel_path}")
print(f"Added {len(libraries) * 6} Mature/Star columns, 6 Max columns (4 values + 2 library names), and removed {len(cols_to_drop)} legacy *_percent columns")
print(f"Enhanced Excel has {len(enhanced_df)} rows and {len(enhanced_df.columns)} columns")

# Print summary of added columns
added_columns = [col for col in enhanced_df.columns if col.endswith('_5p_percent') or col.endswith('_3p_percent')]
print(f"Added columns: {len(added_columns)} total")
print("Sample of added columns:", added_columns[:10])

print("FINISH enhancing Excel file")