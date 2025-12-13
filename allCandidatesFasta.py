import pandas as pd
import sys


def create_all_candidatess_fasta(df):
    """
    Generates FASTA files for mature, star, and hairpin sequences from a DataFrame.
    It creates two versions of the mature and star files: one with the original
    description, and a second version with a '-5p' or '-3p' suffix added to the
    description to indicate the arm of origin, and with 'U' replaced by 'T'.
    It also creates a second version of the hairpin file with 'U' replaced by 'T'.
    """
    # Define output file paths
    fasta_path = output + "all_candidates_mature.fasta"
    fasta_star_path = output + "all_candidates_star.fasta"
    fasta_pre_only_path = output + "all_candidates_hairpin.fasta"

    # --- MODIFICATION START ---
    # Define paths for the new files with 5p/3p suffixes and T-substitution
    fasta_p_path = output + "all_candidates_mature_p.fasta"
    fasta_star_p_path = output + "all_candidates_star_p.fasta"
    fasta_pre_only_T_path = output + "all_candidates_hairpin_T.fasta"

    # Consolidate all file paths to clear them efficiently
    all_files = [
        fasta_path, fasta_star_path, fasta_pre_only_path,
        fasta_p_path, fasta_star_p_path, fasta_pre_only_T_path
    ]

    # Clear all output files before writing
    for f_path in all_files:
        open(f_path, 'w').close()

    # Initialize string buffers for accumulating FASTA content
    fasta_file = ''
    fasta_star_file = ''
    fasta_pre_only_file = ''
    fasta_p_file = ''
    fasta_star_p_file = ''
    fasta_pre_only_T_file = ''
    # --- MODIFICATION END ---

    for index, row in df.iterrows():
        seq5p = row['5pseq']
        mature_side = row['mature']
        seq3p = row['3pseq']
        hairpin = row['hairpinSeq']
        name = row['Description']

        # Process rows where the mature sequence is on the 5p arm
        if not pd.isnull(seq5p) and not pd.isnull(seq3p) and mature_side == '5p':
            # Original files
            fasta_file += f'>{name}\n{seq5p}\n'
            fasta_star_file += f'>{name}\n{seq3p}\n'
            fasta_pre_only_file += f'>{name}\n{hairpin}\n'

            # --- MODIFICATION START ---
            # New files with -5p/-3p suffix and T-substitution
            fasta_p_file += f'>{name}-5p\n{seq5p.replace("U", "T")}\n'
            fasta_star_p_file += f'>{name}-3p\n{seq3p.replace("U", "T")}\n'
            fasta_pre_only_T_file += f'>{name}\n{hairpin.replace("U", "T")}\n'
            # --- MODIFICATION END ---

        # Process rows where the mature sequence is on the 3p arm
        if not pd.isnull(seq3p) and not pd.isnull(seq5p) and mature_side == '3p':
            # Original files
            fasta_file += f'>{name}\n{seq3p}\n'
            fasta_star_file += f'>{name}\n{seq5p}\n'
            fasta_pre_only_file += f'>{name}\n{hairpin}\n'

            # --- MODIFICATION START ---
            # New files with -5p/-3p suffix and T-substitution
            fasta_p_file += f'>{name}-3p\n{seq3p.replace("U", "T")}\n'
            fasta_star_p_file += f'>{name}-5p\n{seq5p.replace("U", "T")}\n'
            fasta_pre_only_T_file += f'>{name}\n{hairpin.replace("U", "T")}\n'
            # --- MODIFICATION END ---

        # Periodically write buffer to files to manage memory usage
        if len(fasta_file) > 100000:
            with open(fasta_path, 'a') as f:
                f.write(fasta_file)
            fasta_file = ''

        if len(fasta_star_file) > 100000:
            with open(fasta_star_path, 'a') as f:
                f.write(fasta_star_file)
            fasta_star_file = ''

        if len(fasta_pre_only_file) > 100000:
            with open(fasta_pre_only_path, 'a') as f:
                f.write(fasta_pre_only_file)
            fasta_pre_only_file = ''

        # --- MODIFICATION START ---
        # Buffer writing for the new files
        if len(fasta_p_file) > 100000:
            with open(fasta_p_path, 'a') as f:
                f.write(fasta_p_file)
            fasta_p_file = ''

        if len(fasta_star_p_file) > 100000:
            with open(fasta_star_p_path, 'a') as f:
                f.write(fasta_star_p_file)
            fasta_star_p_file = ''

        if len(fasta_pre_only_T_file) > 100000:
            with open(fasta_pre_only_T_path, 'a') as f:
                f.write(fasta_pre_only_T_file)
            fasta_pre_only_T_file = ''
        # --- MODIFICATION END ---

    # Write any remaining content from the buffers to the files
    with open(fasta_path, 'a') as f:
        f.write(fasta_file)
    with open(fasta_star_path, 'a') as f:
        f.write(fasta_star_file)
    with open(fasta_pre_only_path, 'a') as f:
        f.write(fasta_pre_only_file)

    # --- MODIFICATION START ---
    # Final write for the new files
    with open(fasta_p_path, 'a') as f:
        f.write(fasta_p_file)
    with open(fasta_star_p_path, 'a') as f:
        f.write(fasta_star_p_file)
    with open(fasta_pre_only_T_path, 'a') as f:
        f.write(fasta_pre_only_T_file)
    # --- MODIFICATION END ---


if __name__ == '__main__':
    print("START")
    species = None
    all_path = None
    sheet_name = "all_candidates"
    new_genome = False
    # Set a default output path
    output = "./"
    output_explicitly_provided = False

    # Parse command-line arguments
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg == '--new-genome':
            new_genome = True
            i += 1
            continue
        if arg == '-s':
            species = sys.argv[i + 1]
        elif arg == '--all':
            all_path = sys.argv[i + 1]
        elif arg == "--sheetname":
            sheet_name = sys.argv[i + 1]
        elif arg == "--output":
            output = sys.argv[i + 1]
            output_explicitly_provided = True
            if not output.endswith('/'):
                output += '/'
        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -s <name>: name of species.\n'
                  f' --all <path>: path to an intersection table excel which contains "all_candidates" sheet.\n'
                  f' --sheetname <name>: name of the sheet in the excel file (default: all_candidates).\n'
                  f' --output <path>: path to the output directory (default: ./).\n'
                  f' --new-genome: use new genome folder structure for output (overrides default ./ if --output not specified).\n'
                  )
            sys.exit()
        i += 2

    # Override output path if --new-genome flag is set and --output was not explicitly provided
    if new_genome and not output_explicitly_provided:
        output = "/groups/vaksler-group/IsanaRNA/Isana_Tzah/RNAcentral/miRNAs/Hofstenia_newGenome/"

    if not all_path:
        print("Error: Missing required argument --all <path>")
        sys.exit(1)

    print(f"Reading sheet: {sheet_name}")
    all_df = pd.read_excel(all_path, sheet_name=sheet_name)

    # Clean and merge description columns
    if species and species.lower() == "elegans":
        all_df['Description'] = all_df[['Description_mirdeep', 'Description_sRNAbench', 'Description_mirbase']].astype(
            str).agg('__'.join, axis=1).str.replace("|", ",", regex=False).str.replace('.', '_', regex=False)
    else:
        all_df['Description'] = all_df[['Description_mirdeep', 'Description_sRNAbench']].astype(str).agg('__'.join,
                                                                                                         axis=1).str.replace(
            ";", "|", regex=False).str.replace("ID=", "", regex=False).str.replace('.', '', regex=False)

    # Remove nan values from description
    all_df['Description'] = all_df['Description'].str.replace('nan, ', '').str.replace(', nan', '')

    create_all_candidatess_fasta(all_df)
    print("FINISH")
