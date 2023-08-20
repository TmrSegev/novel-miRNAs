import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def count_non_nan_columns(row):
    return row.count()

def save_back_to_all(sheet_dict, filepath):
    # Save to original all candidates files
    sheet_dict["(D) Structural Features"].drop(['Species'], axis=1, inplace=True)

    # Merge seeds by species columns
    overlap_columns = sheet_dict["(D) Structural Features"].columns.intersection(seeds_by_species.columns)
    overlap_columns = overlap_columns.drop(["Seed", "Family"])
    sheet_dict["(D) Structural Features"] = pd.merge(sheet_dict["(D) Structural Features"], seeds_by_species.drop(columns=overlap_columns), on=["Seed", "Family"], how="left")

    # Save the changes back to the same sheet
    with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
        writer.book = writer.book  # Needed for openpyxl compatibility
        writer.sheets = {ws.title: ws for ws in writer.book.worksheets}
        # elegans.to_excel(writer, "(D) Structural Features", index=False)
        for sheet_name, df in sheet_dict.items():
            df.to_excel(writer, sheet_name, index=False)

xls = pd.ExcelFile("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx")
elegans_sheet_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}

xls = pd.ExcelFile("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Macrosperma.xlsx")
macrosperma_sheet_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}

xls = pd.ExcelFile("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Sulstoni.xlsx")
sulstoni_sheet_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}

# elegans = elegans_sheet_dict["(D) Structural Features"]
#elegans = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx", sheet_name="(D) Structural Features")
# macrosperma = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Macrosperma.xlsx", sheet_name="(D) Structural Features")
# sulstoni = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Sulstoni.xlsx", sheet_name="(D) Structural Features")

elegans_sheet_dict["(D) Structural Features"]['Species'] = "Elegans"
macrosperma_sheet_dict["(D) Structural Features"]['Species'] = "Macrosperma"
sulstoni_sheet_dict["(D) Structural Features"]['Species'] = "Sulstoni"

all = pd.concat([elegans_sheet_dict["(D) Structural Features"], macrosperma_sheet_dict["(D) Structural Features"], sulstoni_sheet_dict["(D) Structural Features"]])
all.to_excel("all_species_candidates.xlsx", index=False)
all['Description'] = all[['Description_mirdeep', 'Description_sRNAbench']].astype(str).agg(', '.join, axis=1)
seeds_by_species = pd.pivot_table(all, values='Description', index=['Seed', 'Family'], columns='Species', aggfunc='count')
seeds_by_species['Interspecies_seed_count'] = seeds_by_species.apply(count_non_nan_columns, axis=1)
seeds_by_species.fillna(0, inplace=True)
seeds_by_species.to_csv("seeds_by_species.csv")
seeds_by_species = pd.read_csv("seeds_by_species.csv")
known = seeds_by_species[seeds_by_species['Family'] != "UNKNOWN"].copy()
unknown = seeds_by_species[seeds_by_species['Family'] == "UNKNOWN"].copy()
unknown_shared = unknown[unknown['Interspecies_seed_count'] > 1].copy()
unknown_unique = unknown[unknown['Interspecies_seed_count'] == 1].copy()

X_axis = np.arange(3)
# known['Count'].value_counts().plot.bar(X_axis - 0.2, 0.4, label='Known', color="orange")
# unknown['Count'].value_counts().plot.bar(X_axis + 0.2, 0.4, label='Unknown', color="blue")
plt.bar(X_axis - 0.2, known['Interspecies_seed_count'].value_counts(), 0.4, label='Known')
plt.bar(X_axis + 0.2, unknown['Interspecies_seed_count'].value_counts(), 0.4, label='Unknown')
plt.xticks(X_axis, [1, 2, 3])
plt.xlabel("Number of species sharing the seed")
plt.ylabel("Number of Seeds")
plt.title("Known/Unknown Seeds between Species")
plt.legend()
plt.savefig("known_unknown_seeds_between_species.png", dpi=300)
# seeds_by_species.drop("Family", axis=1, inplace=True)
seeds_by_species.to_csv("seeds_by_species.csv", index=False)

unknown_shared_seeds = unknown_shared['Seed'].to_list()
all_unknown_shared = all[all['Seed'].isin(unknown_shared_seeds)]
all_unknown_shared.to_csv("all_unknown_shared.csv", index=False)

def autopct_format(values):
    def my_format(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return 'n={v:d}'.format(v=val)
    return my_format


plt.clf()
unknown_unique['Sum'] = unknown_unique["Elegans"] + unknown_unique["Macrosperma"] + unknown_unique["Sulstoni"]
unknown_unique['Sum'] = unknown_unique['Sum'].astype('int64')
# unknown_unique['Sum'].plot.hist()
# plt.savefig("candidate_counts_of_unknown_seeds_in_one_species.png", dpi=300)
value_counts = unknown_unique['Sum'].value_counts().sort_index()
value_counts.plot.pie(figsize=(9, 9), autopct=autopct_format(value_counts.values), pctdistance=0.9, radius=1.2)
plt.legend()
plt.ylabel("")
plt.title("Candidate counts of unknown seeds that are present in any one species")
plt.savefig("candidate_counts_of_unknown_seeds_in_any_one_species.png", dpi=300)


plt.clf()
value_counts = unknown_unique[unknown_unique["Elegans"] > 0]['Sum'].value_counts().sort_index()
value_counts.plot.pie(figsize=(9, 9), autopct=autopct_format(value_counts.values), pctdistance=0.9, radius=1.2)
plt.legend()
plt.ylabel("")
plt.title("Candidate counts of unknown seeds that are present only in Elegans")
plt.savefig("candidate_counts_of_unknown_seeds_in_elegans.png", dpi=300)

plt.clf()
value_counts = unknown_unique[unknown_unique["Macrosperma"] > 0]['Sum'].value_counts().sort_index()
value_counts.plot.pie(figsize=(9, 9), autopct=autopct_format(value_counts.values), pctdistance=0.9, radius=1.2)
plt.legend()
plt.ylabel("")
plt.title("Candidate counts of unknown seeds that are present only in Macrosperma")
plt.savefig("candidate_counts_of_unknown_seeds_in_macrosperma.png", dpi=300)

plt.clf()
value_counts = unknown_unique[unknown_unique["Sulstoni"] > 0]['Sum'].value_counts().sort_index()
value_counts.plot.pie(figsize=(9, 9), autopct=autopct_format(value_counts.values), pctdistance=0.9, radius=1.2)
plt.legend()
plt.ylabel("")
plt.title("Candidate counts of unknown seeds that are present only in Sulstoni")
plt.savefig("candidate_counts_of_unknown_seeds_in_sulstoni.png", dpi=300)

unknown_unique_two_plus = unknown_unique[unknown_unique['Sum'] > 1]
unknown_unique_two_plus_seeds = unknown_unique_two_plus['Seed'].to_list()
all_unknown_unique_two_plus = all[all['Seed'].isin(unknown_unique_two_plus_seeds)]
all_unknown_unique_two_plus.to_excel("all_unknown_unique_two_plus.xlsx", index=False)

save_back_to_all(elegans_sheet_dict, "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx")
save_back_to_all(macrosperma_sheet_dict, "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Macrosperma.xlsx")
save_back_to_all(sulstoni_sheet_dict, "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Sulstoni.xlsx")
# # Save to original all candidates files
# elegans_sheet_dict["(D) Structural Features"].drop(['Species'], axis=1, inplace=True)
# macrosperma.drop(['Species'], axis=1, inplace=True)
# sulstoni.drop(['Species'], axis=1, inplace=True)
#
# # Merge seeds by species columns
# elegans_sheet_dict["(D) Structural Features"] = pd.merge(elegans_sheet_dict["(D) Structural Features"], seeds_by_species, on=["Seed", "Family"], how="left")
# macrosperma = pd.merge(macrosperma, seeds_by_species, on=["Seed", "Family"], how="left")
# sulstoni = pd.merge(sulstoni, seeds_by_species, on=["Seed", "Family"], how="left")
#
# # Save the changes back to the same sheet
# with pd.ExcelWriter("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx", engine='openpyxl') as writer:
#     writer.book = writer.book  # Needed for openpyxl compatibility
#     writer.sheets = {ws.title: ws for ws in writer.book.worksheets}
#     # elegans.to_excel(writer, "(D) Structural Features", index=False)
#     print(elegans_sheet_dict.items())
#     for sheet_name, df in elegans_sheet_dict.items():
#         print(sheet_name)
#         print(df)
#         df.to_excel(writer, sheet_name, index=False)
# elegans.to_excel("elegans_all_candidates.xlsx", sheet_name='all_candidates')
# macrosperma.to_excel("macrosperma_all_candidates.xlsx", sheet_name='all_candidates')
# sulstoni.to_excel("sulstoni_all_candidates.xlsx", sheet_name='all_candidates')