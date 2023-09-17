import matplotlib.pyplot as plt
import pandas as pd
import sys
from Bio import AlignIO
import os
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Phylo



# def create_graph(row):
#     print(row)
#     plt.clf()
#     row_rpm = row[time]
#     row_rpm.plot.line(style='.-')
#     plot = plt.ylabel("Mature Counts (RPM)")
#     # plots.append(plot)
#     plt.xlabel("Time (Hours)")
#     # plt.title(row['Description'])
#     plt.xticks(ticks=range(0, len(time)), labels=time)
#     plt.savefig('{}_{}.png'.format(row['Seed'], row['Unnamed: 0']), dpi=300)
#
def create_seed_fasta(df):
    fasta_path = "./{}/{}_{}.fasta".format(folder_name, seed, family)
    fasta_file = ''
    open(fasta_path, 'w').close()

    for index, row in df.iterrows():
        species = row['Species']
        hairpin = row['hairpinSeq']
        if not pd.isnull(row['subject_accession']):
            name = row['subject_accession']
        else:
            name = row['Description_mirbase'].split(';')
            name = name[0] + ';' + name[1]
        fasta_file += f'>index={index}|{species}|{name}\n{hairpin}\n'

        if len(fasta_file) > 100000:
            with open(fasta_path, 'a+') as f:
                f.write(fasta_file)
            fasta_file = ''

    with open(fasta_path, 'a+') as f:
        f.write(fasta_file)

def create_folder(seed, family):
    folder_name = "{}_{}".format(seed, family)
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
        print(f"Folder '{folder_name}' created successfully.")
    else:
        print(f"Folder '{folder_name}' already exists.")
    return folder_name

all_path = None
seed = None
libraries_elegans_rpm = ['CE81_m_rpm', 'CE80_m_rpm', 'CE69_m_rpm', 'CE63_m_rpm', 'CE62_m_rpm', 'CE61_m_rpm', 'CE60_m_rpm', 'CE59_m_rpm', 'CE58_m_rpm', 'CE57_m_rpm', 'CE79_m_rpm', 'CE78_m_rpm']
time_elegans = ['4', '8', '12', '16', '20', '24', '28', '32', '36', '40', '44', '48']
libraries_macrosperma_rpm = ['MR8_m_rpm', 'MR7_m_rpm', 'MR6_m_rpm', 'MR5_m_rpm', 'MR4_m_rpm']
time_macrosperma = ['8', '17', '22', '29', '33']
libraries_sulstoni_rpm = ['SR7_m_rpm', 'SR6_m_rpm', 'SR5_m_rpm', 'SR4_m_rpm', 'SR3_m_rpm', 'SR2_m_rpm', 'SR1_m_rpm', 'SR0_m_rpm']
time_sulstoni = ['4', '8', '12', '16', '20', '24', '28', '32']
for i in range(1, len(sys.argv), 2):
    arg = sys.argv[i]
    if arg == '--all':
        all_path = sys.argv[i + 1]
    elif arg == '--seed':
        seed = sys.argv[i + 1]
    elif arg == '--help' or arg == '-h':
        print(f'Manual:\n'
              f' --all <path>: path to an intersection table excel which contains "all_candidates" sheet.\n'
              )
        sys.exit()

# xls = pd.ExcelFile(all_path)
# sheet_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}
all_candidates = pd.read_excel(all_path)
# all_candidates = sheet_dict[0]
# libraries_rpm_all = libraries_elegans_rpm + libraries_macrosperma_rpm + libraries_sulstoni_rpm
# time_all = time_elegans + time_macrosperma + time_sulstoni
# rename_dict = dict(zip(libraries_rpm_all, time_all))
# all_candidates = all_candidates.rename(columns=rename_dict)
seed_candidates = all_candidates[all_candidates['Seed'] == seed]
seed_candidates = seed_candidates.reset_index(drop=True)
family = seed_candidates['Family'].iloc[0]
folder_name = create_folder(seed, family)
create_seed_fasta(seed_candidates)

# plots = []
# seed_candidates.apply(lambda row: create_graph(row), axis=1)
plt.figure(figsize=(22, 14))
plt.suptitle('Seed: {}, Family: {}'.format(seed_candidates['Seed'].iloc[0], seed_candidates['Family'].iloc[0]), fontsize=16)
i = 1
j = ""
seed_candidates.to_excel("./{}/{}_{}_candidates.xlsx".format(folder_name, seed, family), index=False)

for index, row in seed_candidates.iterrows():
    plt.subplot(3, 4, i)
    species = row['Species']
    if species == "Elegans":
        rpm = libraries_elegans_rpm
        time = time_elegans
    elif species == "Macrosperma":
        rpm = libraries_macrosperma_rpm
        time = time_macrosperma
    elif species == "Sulstoni":
        rpm = libraries_sulstoni_rpm
        time = time_sulstoni
    row_rpm = row[rpm]
    rename_dict = dict(zip(rpm, time))
    row_rpm = row_rpm.rename(rename_dict)
    row_rpm.plot.line(style='.-')
    plot = plt.ylabel("Mature Counts (RPM)")
    # plots.append(plot)
    plt.xlabel("Time (Hours)")
    plt.title("index: " + str(index) + "|Species: " + str(species))
    plt.xticks(ticks=range(0, len(time)), labels=time)
    i += 1
    if i == 13:
        i = 1
        if j == "":
            j = 1
        plt.savefig('./{}/{}_{}_{}.png'.format(folder_name, row['Seed'], row['Family'], j), dpi=300)
        plt.clf()
        j += 1

if j == "":
    plt.savefig('./{}/{}_{}.png'.format(folder_name, row['Seed'], row['Family']), dpi=300)
else:
    plt.savefig('./{}/{}_{}_{}.png'.format(folder_name, row['Seed'], row['Family'], j), dpi=300)

file_name = "./{}/{}_{}".format(folder_name, seed, family)
clustalw_exe = r"/sise/home/stome/clustaw/clustalw-2.1/src/clustalw2"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=file_name + ".fasta")
assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
stdout, stderr = clustalw_cline()
alignment = AlignIO.read(file_name+".aln", "clustal")
print(alignment)

calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)
print(distance_matrix)

# ------ Using UPGMA algorithm
constructor = DistanceTreeConstructor()
# Construct the phlyogenetic tree using UPGMA algorithm
UPGMATree = constructor.upgma(distance_matrix)

# Make a better looking tree using the features of matplotlib
fig = plt.figure(figsize=(20, 12), dpi=300) # create figure & set the size
plt.rc('font', size=6)              # fontsize of the leaf and node labels
plt.rc('xtick', labelsize=8)       # fontsize of the tick labels
plt.rc('ytick', labelsize=8)       # fontsize of the tick labels
axes = fig.add_subplot(1, 1, 1)

# drawing the tree
Phylo.draw(UPGMATree, axes=axes)
fig.savefig(file_name + "_UPGMATree.png")

# Further info and find common ancestor of two terminals:
# UPGMATree.common_ancestor('Dog','Masked_Palm_Civet')
# UPGMATree.get_terminals()
# UPGMATree.get_nonterminals()
# UPGMATree.count_terminals()

# ------ Using Neighbour Joining algorithm
# Construct the phlyogenetic tree using NJ algorithm
NJTree = constructor.nj(distance_matrix)
# Draw the phlyogenetic tree using terminal
Phylo.draw_ascii(NJTree)
# Make a better looking tree using the features of matplotlib

fig = plt.figure(figsize=(28, 11), dpi=300) # create figure & set the size
plt.rc('font', size=6)              # fontsize of the leaf and node labels
# matplotlib.rc('font', size=18)              # fontsize of the leaf and node labels
# matplotlib.rc('xtick', labelsize=16)       # fontsize of the tick labels
# matplotlib.rc('ytick', labelsize=16)       # fontsize of the tick labels

axes = fig.add_subplot(1, 1, 1)
Phylo.draw(NJTree, axes=axes)

fig.savefig(file_name + "_NJTree.png")

