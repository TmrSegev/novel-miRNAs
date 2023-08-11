import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def count_non_nan_columns(row):
    return row.count()

def save_back_to_all(species, filepath):


xls = pd.ExcelFile("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans_backup.xlsx")
elegans_sheet_dict = {sheet_name: xls.parse(sheet_name) for sheet_name in xls.sheet_names}
print(elegans_sheet_dict)
# elegans = elegans_sheet_dict["(D) Structural Features"]
#elegans = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx", sheet_name="(D) Structural Features")
macrosperma = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Macrosperma.xlsx", sheet_name="(D) Structural Features")
sulstoni = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Sulstoni.xlsx", sheet_name="(D) Structural Features")

elegans_sheet_dict["(D) Structural Features"]['Species'] = "Elegans"
macrosperma['Species'] = "Macrosperma"
sulstoni['Species'] = "Sulstoni"

all = pd.concat([elegans_sheet_dict["(D) Structural Features"], macrosperma, sulstoni])
all['Description'] = all[['Description_mirdeep', 'Description_sRNAbench']].astype(str).agg(', '.join, axis=1)
seeds_by_species = pd.pivot_table(all, values='Description', index=['Seed', 'Family'], columns='Species', aggfunc='count')
seeds_by_species['Count'] = seeds_by_species.apply(count_non_nan_columns, axis=1)
seeds_by_species.fillna(0, inplace=True)
seeds_by_species.to_csv("seeds_by_species.csv")
seeds_by_species = pd.read_csv("seeds_by_species.csv")
known = seeds_by_species[seeds_by_species['Family'] != "UNKNOWN"].copy()
unknown = seeds_by_species[seeds_by_species['Family'] == "UNKNOWN"].copy()
unknown_shared = unknown[unknown['Count'] > 1].copy()

X_axis = np.arange(3)
# known['Count'].value_counts().plot.bar(X_axis - 0.2, 0.4, label='Known', color="orange")
# unknown['Count'].value_counts().plot.bar(X_axis + 0.2, 0.4, label='Unknown', color="blue")
plt.bar(X_axis - 0.2, known['Count'].value_counts(), 0.4, label='Known')
plt.bar(X_axis + 0.2, unknown['Count'].value_counts(), 0.4, label='Unknown')
plt.xticks(X_axis, [1, 2, 3])
plt.xlabel("Number of species sharing the seed")
plt.ylabel("Number of Seeds")
plt.title("Known/Unknown Seeds between Species")
plt.legend()
plt.savefig("known_unknown_seeds_between_species.png", dpi=300)
# seeds_by_species.drop("Family", axis=1, inplace=True)
seeds_by_species.to_csv("seeds_by_species.csv", index=False)

# Save to original all candidates files
elegans_sheet_dict["(D) Structural Features"].drop(['Species'], axis=1, inplace=True)
macrosperma.drop(['Species'], axis=1, inplace=True)
sulstoni.drop(['Species'], axis=1, inplace=True)

# Merge seeds by species columns
elegans_sheet_dict["(D) Structural Features"] = pd.merge(elegans_sheet_dict["(D) Structural Features"], seeds_by_species, on=["Seed", "Family"], how="left")
macrosperma = pd.merge(macrosperma, seeds_by_species, on=["Seed", "Family"], how="left")
sulstoni = pd.merge(sulstoni, seeds_by_species, on=["Seed", "Family"], how="left")

# Save the changes back to the same sheet
with pd.ExcelWriter("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx", engine='openpyxl') as writer:
    writer.book = writer.book  # Needed for openpyxl compatibility
    writer.sheets = {ws.title: ws for ws in writer.book.worksheets}
    # elegans.to_excel(writer, "(D) Structural Features", index=False)
    print(elegans_sheet_dict.items())
    for sheet_name, df in elegans_sheet_dict.items():
        print(sheet_name)
        print(df)
        df.to_excel(writer, sheet_name, index=False)
    # Writer.save()
# elegans.to_excel("elegans_all_candidates.xlsx", sheet_name='all_candidates')
# macrosperma.to_excel("macrosperma_all_candidates.xlsx", sheet_name='all_candidates')
# sulstoni.to_excel("sulstoni_all_candidates.xlsx", sheet_name='all_candidates')