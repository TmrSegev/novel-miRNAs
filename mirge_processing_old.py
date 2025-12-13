import os
import pandas as pd
import re
import argparse
from pathlib import Path

# --- Argument Parsing ---
# Set up a parser to read command-line arguments
parser = argparse.ArgumentParser(
    description=""
)
# Add a required positional argument for the base directory path
parser.add_argument("--dir", help="")
parser.add_argument("--species", help="")
# Parse the arguments provided by the user
args = parser.parse_args()
dir = args.dir
SPECIES = args.species

# Configure the input Excel file path and sheet name
input_excel = Path(f"/groups/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_{SPECIES}.xlsx")
# Note: Corrected a potential typo in the sheet name from "(A) Unfiltered)"
sheet_name = "(A) Unfiltered" if SPECIES == "Hofstenia" else "(D) Structural Features"

print(f"Loading maturity info from: {input_excel} | Sheet: {sheet_name}")
# Load the maturity data from the specified sheet in the Excel file
maturity_df = pd.read_excel(input_excel, sheet_name=sheet_name)

# Create a mapping dictionary from precursor miRNA name ('Description') to its mature arm ('mature')
# This will be used to determine if a miRNA is the mature or star sequence.
maturity_df = maturity_df[['Description', 'mature']].dropna()
maturity_map = pd.Series(maturity_df['mature'].values, index=maturity_df['Description'].str.strip()).to_dict()

# Define the base directory containing all libraries
# base_dir = '/groups/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/results/'
base_dir = dir + '/results_m18/'

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
        total_hits_filtered = group[~group['Variant'].str.contains('iso_add5p|iso_add3p', na=False)]['Hits'].sum()
        ref_hits = group[group['Type'] == 'ref_miRNA']['Hits'].sum()

        higher_isomir = 'NA'
        higher_isomir_seq = 'NA'
        higher_isomir_length = 'NA'
        higher_isomir_count = 'NA'

        isomir_group = group[group['Type'] == 'isomiR']
        if not isomir_group.empty:
            max_isomir_row = isomir_group.loc[isomir_group['Hits'].idxmax()]
            if ref_hits == 0 or max_isomir_row['Hits'] > ref_hits:
                higher_isomir = max_isomir_row['Variant']
                higher_isomir_seq = max_isomir_row['Sequence']
                higher_isomir_length = max_isomir_row['Length']
                higher_isomir_count = max_isomir_row['Hits']

        iso_5p_hits = group[group['Variant'].str.contains(r'iso_5p:', na=False) &
                            ~group['Variant'].str.contains('iso_add5p|iso_add3p', na=False)]['Hits'].sum()
        iso_3p_hits = group[group['Variant'].str.contains(r'iso_3p', na=False) &
                            ~group['Variant'].str.contains('iso_add5p|iso_add3p', na=False)]['Hits'].sum()
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

    # # Save the individual library report
    # output_dir = "/groups/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/processing_results/"
    # output_filename = os.path.join(output_dir, f'miRNA_isomiR_report_{library_name}.csv')
    # report_df.to_csv(output_filename, index=False)

    # # Pairing 5p and 3p miRNAs
    # def extract_pre_mirna_name(mirna_name):
    #     return mirna_name.rsplit('-', 1)[0]
    #
    # paired_data = []
    # feature_columns = [col for col in report_df.columns if col not in ['miRNA', 'Library']]
    #
    # for pre_mirna, group in report_df.groupby(report_df['miRNA'].apply(extract_pre_mirna_name)):
    #     group = group.set_index('miRNA')
    #
    #     if any('-5p' in mir for mir in group.index) and any('-3p' in mir for mir in group.index):
    #         miRNA_5p = group.filter(like='-5p', axis=0)
    #         miRNA_3p = group.filter(like='-3p', axis=0)
    #
    #         if not miRNA_3p.empty and not miRNA_5p.empty:
    #             row = [pre_mirna]
    #             row.extend(miRNA_3p.iloc[0][feature_columns].tolist())
    #             row.extend(miRNA_5p.iloc[0][feature_columns].tolist())
    #             row.append(miRNA_3p.iloc[0]['Library'])
    #             paired_data.append(row)
    #     else:
    #         no_3p_or_5p_mirs_overall += 1
    #         no_3p_or_5p_mirs_library += 1
    #
    # paired_columns = ['pre-miRNA'] + [f'3p_{col}' for col in feature_columns] + [f'5p_{col}' for col in feature_columns] + ['Library']
    # paired_df = pd.DataFrame(paired_data, columns=paired_columns)
    #
    # # Store all paired data
    # all_paired_data.append(paired_df)
    #
    # print(f"MIRs removed for not having 3p or 5p in library {library_name}: {no_3p_or_5p_mirs_library}")

# Concatenate all paired data into a final dataframe
# final_paired_df = pd.concat(all_paired_data, ignore_index=True)
final_df = pd.concat(all_data, ignore_index=True)

# Sort by pre-miRNA and Library
# final_paired_df.sort_values(by=['pre-miRNA', 'Library'], inplace=True)
# final_df.sort_values(by=['pre-miRNA', 'Library'], inplace=True)
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


# Apply the function to each row to create the new 'mature' column
final_df['mature'] = final_df.apply(lambda row: get_maturity_status(row, maturity_map), axis=1)

# Reorder columns to place 'mature' right after 'miRNA' for better readability
cols = final_df.columns.tolist()
# Pop the 'mature' column from the end and insert it at the second position
cols.insert(1, cols.pop(cols.index('mature')))
final_df = final_df[cols]
# --- END NEW SECTION ---

# Save the final concatenated and sorted paired dataframe
# output_dir = "/groups/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/miRge_output/processing_results/"
output_dir = dir + "/processing_results_m18/"
os.makedirs(output_dir, exist_ok=True)
# final_output_filename = os.path.join(output_dir, 'final_miRNA_isomiR_paired_report.csv')
# final_paired_df.to_csv(final_output_filename, index=False)
final_output_filename = os.path.join(output_dir, 'final_miRNA_isomiR_report_m18.csv')
final_df.to_csv(final_output_filename, index=False)

# print("Paired report saved for all libraries.")
# print(f"MIRs removed for not having 3p or 5p across all libraries: {no_3p_or_5p_mirs_overall}")
print("FINISH")
