import matplotlib.pyplot as plt
import pandas as pd
import sys


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
    fasta_path = "{}.fasta".format(seed)
    fasta_file = ''
    open(fasta_path, 'w').close()

    for index, row in df.iterrows():

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


species = None
all_path = None
seed = None
libraries = None
time=None
# libraries_elegans_rpm = ['CE81_m_rpm', 'CE80_m_rpm', 'CE69_m_rpm', 'CE63_m_rpm', 'CE62_m_rpm', 'CE61_m_rpm', 'CE60_m_rpm', 'CE59_m_rpm', 'CE58_m_rpm', 'CE57_m_rpm', 'CE79_m_rpm', 'CE78_m_rpm']
# time_elegans = [4,8,12,16,20,24,28,32,36,40,44,48]
# libraries_macrosperma_rpm = ['MR8_m_rpm', 'MR7_m_rpm', 'MR6_m_rpm', 'MR5_m_rpm', 'MR4_m_rpm']
# time_macrosperma = [8, 17, 22, 29, 33]
# libraries_sulstoni_rpm = ['SR7_m_rpm', 'SR6_m_rpm', 'SR5_m_rpm', 'SR4_m_rpm', 'SR3_m_rpm', 'SR2_m_rpm', 'SR1_m_rpm', 'SR0_m_rpm']
# time_sulstoni = [4, 8, 12, 16, 20, 24, 28, 32]
for i in range(1, len(sys.argv), 2):
    arg = sys.argv[i]
    if arg == '-s':
        species = sys.argv[i + 1]
    elif arg == '--all':
        all_path = sys.argv[i + 1]
    elif arg == '--seed':
        seed = sys.argv[i + 1]
    elif arg == '--libraries':
        libraries = sys.argv[i + 1].split(',')
    elif arg == '--time':
        time = sys.argv[i + 1].split(',')
    elif arg == '--help' or arg == '-h':
        print(f'Manual:\n'
              f' -s <name>: name of species.\n'
              f' --all <path>: path to an intersection table excel which contains "all_candidates" sheet.\n'
              )
        sys.exit()
xls = pd.ExcelFile(all_path)
sheet_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}
# all = pd.read_excel(all_path, sheet_name="(D) Structural Features")
if species == "all":
    all_candidates = sheet_dict[0]
else:
    all_candidates = sheet_dict["(D) Structural Features"]
libraries_rpm = [library + "_m_rpm" for library in libraries]
rename_dict = dict(zip(libraries_rpm, time))
all_candidates = all_candidates.rename(columns=rename_dict)

seed_candidates = all_candidates[all_candidates['Seed'] == seed]
seed_candidates = seed_candidates.reset_index(drop=True)
create_seed_fasta(seed_candidates)

# plots = []
# seed_candidates.apply(lambda row: create_graph(row), axis=1)
plt.figure(figsize=(12, 14))
plt.suptitle('Seed: {}, Family: {}'.format(seed_candidates['Seed'].iloc[0], seed_candidates['Family'].iloc[0]), fontsize=16)
i = 1
for index, row in seed_candidates.iterrows():
    plt.subplot(3, 2, i)
    row_rpm = row[time]
    row_rpm.plot.line(style='.-')
    plot = plt.ylabel("Mature Counts (RPM)")
    # plots.append(plot)
    plt.xlabel("Time (Hours)")
    plt.title("index: " + str(index))
    plt.xticks(ticks=range(0, len(time)), labels=time)
    i += 1
    if i == 7:
        i = 1
        if j == "":
            j = 1
        plt.savefig('{}_{}{}.png'.format(row['Seed'], row['Family'], j), dpi=300)
        plt.clf()
        j += 1

plt.savefig('{}_{}{}.png'.format(row['Seed'], row['Family'], j), dpi=300)

# DATAQUEST WAY
# plt.figure(figsize=(10, 12))
# for i in range(1, 7):
#     plt.subplot(3, 2, i)
#     plt.plot(plots[i-1])
#
# plt.show()



