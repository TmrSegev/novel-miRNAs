import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
from scipy.stats import ttest_ind

# ---Creating a families by type pivot table


def families_by_type(df):
    families_by_type = pd.pivot_table(df, values='Description', index='Family', columns='Type', aggfunc='count')
    families_by_type.fillna(0, inplace=True)
    families_by_type = families_by_type.drop("UNKNOWN", axis=0)
    families_by_type.sort_values('Family', inplace=True)
    families_by_type.plot(kind='bar', figsize=(14, 10), stacked=True, title="{} counts of known families in each type".format(species))
    plt.ylabel("Counts")
    # plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
    # plt.show()
    plt.savefig("{}_family_counts_per_type.png".format(species))


# ---Creating unknown seeds by type pivot table


def unknown_families_by_type(df):
    filtered_df = df[df['Family'] == "UNKNOWN"]
    pivot = pd.pivot_table(filtered_df, values='Description', index='Seed', columns='Type', aggfunc='count')
    pivot.fillna(0, inplace=True)
    pivot['sum'] = pivot.sum(axis=1)
    pivot_groups = pivot[pivot['sum'] > 1]
    pivot_solos = pivot[pivot['sum'] == 1]

    pivot_groups = pivot_groups.drop('sum', axis=1)
    pivot_groups.sort_values('Seed', inplace=True)
    pivot_groups.plot(kind='bar', figsize=(14, 10), stacked=True,
                          title="{} counts of unknown families in each type".format(species))
    plt.ylabel("Counts")
    # plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
    # plt.show()
    plt.savefig("{}_unknown_family_counts_per_type.png".format(species))
    plt.clf()

    pivot_solos = pivot_solos.drop('sum', axis=1)
    pivot_solos.sum().plot(kind='pie', y='Type', autopct='%1.0f%%',
                           title="{} unique seed candidates by type".format(species))
    plt.ylabel("Counts")
    # plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
    # plt.show()
    plt.savefig("{}_unique_seed_candidates_by_type.png".format(species))


def boxplot_by_type(df):
    plt.clf()
    types = list(df['Type'].unique())
    types.sort()
    columns = []
    for type in types:
        col = 'Type {} (n={})'.format(type, len(df.loc[df['Type'] == type]))
        df[col] = np.log10(df.loc[df['Type'] == type, 'mean_m_rpm'])
        # df[col] = df.loc[df['Type'] == type, 'mean_m_rpm']
        columns.append(col)
    df.boxplot(column=columns)
    plt.ylabel("log10_mean_m_rpm")
    plt.yticks(ticks= np.arange(0, 6), labels=10 ** np.arange(0, 6))
    plt.title("{} boxplots by type".format(species))
    plt.xticks(rotation=45)
    plt.tight_layout()
    print(enumerate(df[columns]))
    for i, d in enumerate(df[columns]):
        y = df[columns][d]
        x = np.random.normal(i + 1, 0.04, len(y))
        plt.scatter(x, y)
    plt.savefig("{}_boxplots_by_type.png".format(species))


def boxplot_known_unknown(df):
    plt.clf()
    types = list(df['Type'].unique())
    types.sort()
    columns = []
    for type in types:
        col_known = 'Type {} known (n={})'.format(type, len(df.loc[(df['Type'] == type) & (df['Family'] != 'UNKNOWN')]))
        col_unknown = 'Type {} unknown (n={})'.format(type, len(df.loc[(df['Type'] == type) & (df['Family'] == 'UNKNOWN')]))
        df[col_known] = df.loc[(df['Type'] == type) & (df['Family'] != 'UNKNOWN'), 'mean_m_rpm']
        df[col_unknown] = df.loc[(df['Type'] == type) & (df['Family'] == 'UNKNOWN'), 'mean_m_rpm']
        df[col_known] = np.log10(df[col_known])
        df[col_unknown] = np.log10(df[col_unknown])
        columns.append(col_known)
        columns.append(col_unknown)
    df.boxplot(column=columns)
    plt.ylabel("log10_mean_m_rpm")
    plt.yticks(ticks=np.arange(0, 6), labels=10 ** np.arange(0, 6))
    plt.title("{} boxplots known/unknown by type".format(species))
    plt.xticks(rotation=45)
    plt.tight_layout()
    for i, d in enumerate(df[columns]):
        y = df[columns][d]
        x = np.random.normal(i + 1, 0.04, len(y))
        plt.scatter(x, y)
    plt.savefig("{}_boxplots_known_unknown.png".format(species))


#def wilcoxon_test(df):
    #print(np.log10(df.loc[df['Type'] == 1, 'mean_m_rpm']).describe())
   # print(np.log10(df.loc[df['Type'] == 1, 'mean_m_rpm']).describe())
  #  stats.wilcoxon(np.log10(df.loc[df['Type'] == 1, 'mean_m_rpm']), np.log10(df.loc[df['Type'] == 2, 'mean_m_rpm']))

def mann_whitney(df):
    types = list(df['Type'].unique())
    types.sort()
    with open('{}_mann_whitney.txt'.format(species), 'w') as file:
        for i in range(1, len(types) + 1):
            for j in range(i + 1, len(types) + 1):
                disti = np.log10(df.loc[df['Type'] == i, 'mean_m_rpm'])
                distj = np.log10(df.loc[df['Type'] == j, 'mean_m_rpm'])
                res = stats.mannwhitneyu(disti, distj, alternative="two-sided")
                file.write("Type " + str(i) + " vs Type " + str(j) + ":    " + str(res) + '\n')

def t_test(df):
    types = list(df['Type'].unique())
    types.sort()
    for i in range(1, len(types)+1):
        for j in range(i+1, len(types)+1):
            disti = np.log10(df.loc[df['Type'] == i, 'mean_m_rpm'])
            distj = np.log10(df.loc[df['Type'] == j, 'mean_m_rpm'])
            stat, p = ttest_ind(disti, distj)
            print("Type", i, "vs Type", j, ":     stat:", stat, ", p:", p)


def normal_dist(df):
    plt.clf()
    types = list(df['Type'].unique())
    types.sort()
    figure, axis = plt.subplots(4, 2)
    for i, ax in enumerate(figure.axes):
        # ax.set_ylabel(str(i))
        ax.hist(df.loc[df['Type'] == type, 'mean_m_rpm'])
        ax.set_title("Type {}".format(type))
    plt.savefig("{}_normal_dist.png".format(species))


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
    families_by_type(all)
    unknown_families_by_type(all)
    boxplot_by_type(all)
    boxplot_known_unknown(all)
    # normal_dist(all)
    # wilcoxon_test(all)
    # t_test(all)
    mann_whitney(all)
