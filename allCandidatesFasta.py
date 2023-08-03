import pandas as pd
import sys

def create_all_candidatess_fasta(df):
    fasta_path = "all_candidates_mature.fasta"
    fasta_star_path = "all_candidates_star.fasta"
    fasta_pre_only_path = "all_candidates_hairpin.fasta"
    fasta_file = ''
    fasta_star_file = ''
    fasta_pre_only_file = ''
    open(fasta_path, 'w').close()
    open(fasta_star_path, 'w').close()
    open(fasta_pre_only_path, 'w').close()

    for index, row in df.iterrows():
        seq5p = row['5pseq']
        mature_seq = row['mature']
        seq3p = row['3pseq']
        hairpin = row['hairpinSeq']
        name = row['Description']
        if not pd.isnull(seq5p) and mature_seq == '5p':
            fasta_file += f'>{name}\n{seq5p}\n'
            fasta_star_file += f'>{name}\n{seq3p}\n'
            fasta_pre_only_file += f'>{name}\n{hairpin}\n'
        if not pd.isnull(seq3p) and mature_seq == '3p':
            fasta_file += f'>{name}\n{seq3p}\n'
            fasta_star_file += f'>{name}\n{seq5p}\n'
            fasta_pre_only_file += f'>{name}\n{hairpin}\n'

        if len(fasta_file) > 100000:
            with open(fasta_path, 'a+') as f:
                f.write(fasta_file)
            fasta_file = ''

        if len(fasta_star_file) > 100000:
            with open(fasta_star_path, 'a+') as f:
                f.write(fasta_star_file)
            fasta_star_file = ''

        if len(fasta_pre_only_file) > 100000:
            with open(fasta_pre_only_path, 'a+') as f:
                f.write(fasta_pre_only_file)
            fasta_pre_only_file = ''

    with open(fasta_path, 'a+') as f:
        f.write(fasta_file)
    with open(fasta_star_path, 'a+') as f:
        f.write(fasta_star_file)
    with open(fasta_pre_only_path, 'a+') as f:
        f.write(fasta_pre_only_file)


if __name__ == '__main__':
    species = None
    all_path = None
    for i in range(1, len(sys.argv), 2):
        arg = sys.argv[i]
        if arg == '-s':
            species = sys.argv[i + 1]
        elif arg == '--all':
            all_path = sys.argv[i + 1]
        elif arg == '--help' or arg == '-h':
            print(f'Manual:\n'
                  f' -s <name>: name of species.\n'
                  f' --all <path>: path to an intersection table excel which contains "all_candidates" sheet.\n'
                  )
            sys.exit()
    all = pd.read_excel(all_path, sheet_name="all_candidates")
    if species == "elegans" or species == "Elegans":
        all['Description'] = all[['Description_mirdeep', 'Description_sRNAbench', 'Description_mirbase']].astype(str).agg(', '.join, axis=1)
    else:
        all['Description'] = all[['Description_mirdeep', 'Description_sRNAbench']].astype(str).agg(', '.join, axis=1)
    # no_novel451 = all[~all['Description'].str.contains("novel451")].copy()
    create_all_candidatess_fasta(all)

