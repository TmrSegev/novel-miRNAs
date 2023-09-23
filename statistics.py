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
    plt.savefig("./figures/{}_family_counts_per_type.png".format(species), dpi=300)


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
    plt.savefig("./figures/{}_unknown_family_counts_per_type.png".format(species), dpi=300)
    plt.clf()

    pivot_solos = pivot_solos.drop('sum', axis=1)
    pivot_solos.sum().plot(kind='pie', y='Type', autopct='%1.0f%%',
                           title="{} unique seed candidates by type".format(species))
    plt.ylabel("Counts")
    # plt.legend(["1: Both", "2: miRDeep only", "3: sRNAbench only"])
    # plt.show()
    plt.savefig("./figures/{}_unique_seed_candidates_by_type.png".format(species), dpi=300)


def boxplot_by_type(df, name):
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
    plt.xticks(rotation=90)
    #plt.tight_layout()
    for i, d in enumerate(df[columns]):
        y = df[columns][d]
        x = np.random.normal(i + 1, 0.04, len(y))
        plt.scatter(x, y)
    plt.savefig("./figures/{}_boxplots_by_type_{}.png".format(species, name), dpi=300, bbox_inches='tight')

    # Create boxplot figure for elegans seperating miRdeep, sRNAbench or both:
    if species == "Elegans" or species == "elegans":
        plt.clf()
        types = [[1, 5, 9], [2, 6, 10], [3, 7, 11]]
        titles = ["Both", "miRdeep only", "sRNAbench only"]
        columns = []
        for i in range(0, len(types)):
            # col = '{} (n={})'.format(titles[i], len(df.loc[df['Type'] in types[i]]))
            # df[col] = np.log10(df.loc[df['Type'] in types[i], 'mean_m_rpm'])
            # # df[col] = df.loc[df['Type'] == type, 'mean_m_rpm']
            # columns.append(col)
            mask = df['Type'].isin(types[i])
            col = '{} (n={})'.format(titles[i], mask.sum())
            df[col] = np.log10(df.loc[mask, 'mean_m_rpm'])
            columns.append(col)
        df.boxplot(column=columns)
        plt.ylabel("log10_mean_m_rpm")
        plt.yticks(ticks=np.arange(0, 6), labels=10 ** np.arange(0, 6))
        plt.title("{} boxplots by miRdeep / sRNAbench".format(species))
        plt.xticks(rotation=90)
        # plt.tight_layout()
        for i, d in enumerate(df[columns]):
            y = df[columns][d]
            x = np.random.normal(i + 1, 0.04, len(y))
            plt.scatter(x, y)
        plt.savefig("./figures/{}_boxplots_by_algorithm_{}.png".format(species, name), dpi=300, bbox_inches='tight')


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
    plt.savefig("./figures/{}_boxplots_known_unknown.png".format(species), dpi=300)


#def wilcoxon_test(df):
    #print(np.log10(df.loc[df['Type'] == 1, 'mean_m_rpm']).describe())
   # print(np.log10(df.loc[df['Type'] == 1, 'mean_m_rpm']).describe())
  #  stats.wilcoxon(np.log10(df.loc[df['Type'] == 1, 'mean_m_rpm']), np.log10(df.loc[df['Type'] == 2, 'mean_m_rpm']))

def mann_whitney(df):
    types = list(df['Type'].unique())
    types.sort()
    with open('./figures/{}_mann_whitney.txt'.format(species), 'w') as file:
        for i in range(1, len(types) + 1):
            for j in range(i + 1, len(types) + 1):
                disti = np.log10(df.loc[df['Type'] == i, 'mean_m_rpm'])
                distj = np.log10(df.loc[df['Type'] == j, 'mean_m_rpm'])
                res_less = stats.mannwhitneyu(disti, distj, alternative="less")
                res_great = stats.mannwhitneyu(disti, distj, alternative="greater")
                min_p = min(res_less.pvalue, res_great.pvalue)
                file.write("Type " + str(i) + " vs Type " + str(j) + ":    " + str(min_p) + '\t' + ("significant" if min_p <= 0.05 else "non-significant") + '\n')

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
    plt.savefig("./figures/{}_normal_dist.png".format(species), dpi=300)

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


def save_back_to_all(sheet_dict, filepath):
    # Save the changes back to the same sheet
    with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
        writer.book = writer.book  # Needed for openpyxl compatibility
        writer.sheets = {ws.title: ws for ws in writer.book.worksheets}
        for sheet_name, df in sheet_dict.items():
            df.to_excel(writer, sheet_name, index=False)

def clusters():
    for key in sheet_dict.keys():
        # Filtering by coordinates
        df = sheet_dict[key].copy()
        df = df.sort_values(['Chr', 'Strand', 'Start', 'End'])
        df['Cluster_index'] = np.zeros(len(df))
        df['Cluster_size'] = np.zeros(len(df))
        clusters_df = pd.DataFrame(columns=df.columns)

        if "blast" not in key:
            cluster_index = 1
            cluster_file = ""

            # Calculate differences within each group
            df['differences'] = df.groupby(['Chr', 'Strand'])['Start'].diff()
            df['differences'] = df['differences'].fillna(0)
            new_df = pd.DataFrame(columns=df.columns)
            # df.to_excel('{}_diff.xlsx'.format(species), index=True)

            # Calculate clusters
            for index, row in df.iterrows():
                if index in df.index:
                    df['distance'] = df['Start'] - row['Start']
                    cluster = df[df['distance'].abs() <= 10000]
                    cluster = cluster[cluster['Chr'] == row['Chr']]
                    cluster = cluster[cluster['Strand'] == row['Strand']]
                    # df.loc[index, 'cluster'] = len(cluster)
                    if len(cluster) == 1:
                        # no_cluster = no_cluster.append(row)
                        cluster['Cluster_index'] = -1
                        cluster['Cluster_size'] = len(cluster)
                        new_df = new_df.append(cluster)
                        df = df.drop(index)
                    else:
                        cluster['Cluster_index'] = cluster_index
                        cluster['Cluster_size'] = len(cluster)
                        cluster_file += "Cluster index:" + str(cluster_index) + "\nSeed families:" + str(cluster["Family"].value_counts().to_dict()) + "\nSeeds:" + str(cluster["Seed"].to_list()) + "\n\n"
                        clusters_df = clusters_df.append(cluster)
                        new_df = new_df.append(cluster)
                        df = df.drop(cluster.index)
                        cluster_index += 1
            # print(df['cluster'].value_counts().sort_index(ascending=False))
            df = df.drop(["distance"], axis=1)
            clusters_df = clusters_df.drop(["distance"], axis=1)
            new_df = new_df.drop(["distance"], axis=1)
            # filtered_input.append(df)
            with open('{}_clusters_info.txt'.format(species), 'w') as file:
                file.write(cluster_file)
            # clusters_df = clusters_df.drop(['Description', 'differences'], axis=1)
            # new_df = new_df.drop(['Description', 'differences'], axis=1)
            new_df.sort_index(inplace=True)
            new_df = new_df.sort_values(['Chr', 'Strand', 'Start', 'End'])
            sheet_dict[key] = new_df
            clusters_df.to_csv('{}_clusters.csv'.format(species), sep='\t', index=False)

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
    xls = pd.ExcelFile(all_path)
    sheet_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}
    # all = pd.read_excel(all_path, sheet_name="(D) Structural Features")
    if species == "Hofstenia":
        all = sheet_dict["(A) Unfiltered)"]
    else:
        all = sheet_dict["(D) Structural Features"]
    if species == "elegans" or species == "Elegans":
        all['Description'] = all[['Description_mirdeep', 'Description_sRNAbench', 'Description_mirbase']].astype(str).agg('. '.join, axis=1)
    else:
        all['Description'] = all[['Description_mirdeep', 'Description_sRNAbench']].astype(str).agg('. '.join, axis=1)
    families_by_type(all.copy())
    unknown_families_by_type(all.copy())
    boxplot_by_type(all.copy(), "all")
    # no_novel451 = all[~all['Description'].str.contains("novel451")].copy()
    # boxplot_by_type(no_novel451, "no_novel451")
    boxplot_known_unknown(all.copy())
    # normal_dist(all)
    # wilcoxon_test(all)
    # t_test(all)
    mann_whitney(all.copy())
    clusters()
    all = all.drop(['Description'], axis=1)
    save_back_to_all(sheet_dict, all_path)
    # create_all_candidatess_fasta(all)
    #pd.to_csv("all_candidates_{}.csv".format(species), sep='\t', index=False)
