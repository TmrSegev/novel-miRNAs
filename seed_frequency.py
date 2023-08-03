import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def count_non_nan_columns(row):
    return row.count()

elegans = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Elegans.xlsx", sheet_name="(D) Structural Features")
macrosperma = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Macrosperma.xlsx", sheet_name="(D) Structural Features")
sulstoni = pd.read_excel("/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Ziv_Features/all_remaining_after_ziv_Sulstoni.xlsx", sheet_name="(D) Structural Features")

elegans['Species'] = "Elegans"
macrosperma['Species'] = "Macrosperma"
sulstoni['Species'] = "Sulstoni"

all = pd.concat([elegans, macrosperma, sulstoni])
all['Description'] = all[['Description_mirdeep', 'Description_sRNAbench']].astype(str).agg(', '.join, axis=1)
seeds_by_species = pd.pivot_table(all, values='Description', index=['Seed', 'Family'], columns='Species', aggfunc='count')
seeds_by_species['Count'] = seeds_by_species.apply(count_non_nan_columns, axis=1)
seeds_by_species.fillna(0, inplace=True)
print(seeds_by_species)
seeds_by_species.to_csv("seeds_by_species.csv")
seeds_by_species = pd.read_csv("seeds_by_species.csv")
print(seeds_by_species)
known = seeds_by_species[seeds_by_species['Family'] != "UNKNOWN"].copy()
unknown = seeds_by_species[seeds_by_species['Family'] == "UNKNOWN"].copy()
print(known['Count'].value_counts())
print(unknown['Count'].value_counts())

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
seeds_by_species.drop("Family", axis=1, inplace=True)
seeds_by_species.to_csv("seeds_by_species.csv")
