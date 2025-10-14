print("START")

import pandas as pd
import argparse
import re

print("IMPORTED")

# Argument parser to toggle between sRNAbench and miRDeep
parser = argparse.ArgumentParser()
parser.add_argument("--tool", choices=["sRNAbench", "miRDeep"], required=True,
                    help="Specify the tool: sRNAbench or miRDeep")
args = parser.parse_args()

tool_name = args.tool
alt_tool_name = "miRDeep" if tool_name == "sRNAbench" else "sRNAbench"

# Load the data
base_path = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/scripts/"
output_path = "/sise/vaksler-group/IsanaRNA/Isana_Tzah/Charles_seq/Hofstenia/good_candidates/"
if tool_name == "miRDeep":
    file1 = f"{base_path}debugging_Hofstenia_miRDeep_1.csv"
    file2 = f"{base_path}debugging_Hofstenia_miRDeep_2.csv"
    df1 = pd.read_csv(file1, sep='\t')
    df2 = pd.read_csv(file2, sep='\t', names=df1.columns, skiprows=1)
    df = pd.concat([df1, df2], ignore_index=True)
    df.to_csv(f'{base_path}debugging_Hofstenia_miRDeep.csv', sep='\t', index=False)
else:
    input_file = f"{base_path}debugging_Hofstenia_{tool_name}.csv"
    df = pd.read_csv(input_file, sep='\t')
    df = df[~df["origin"].str.contains("novel451", na=False)]  # Remove novel451

print("Columns in dataset:", df.columns)

# Adjust column names based on tool
dev_condition_col = "Library"
precursor_col = "hairpinSeq" if tool_name == "sRNAbench" else "consensus precursor sequence"
mature_col = "3pseq" if tool_name == "sRNAbench" else "consensus mature sequence"
star_col = "5pseq" if tool_name == "sRNAbench" else "consensus star sequence"
read_count_col = "totalRC" if tool_name == "sRNAbench" else "total read count"

if tool_name == "miRDeep":
    df[["scaffold", "coordinates", "strand"]] = df["precursor coordinate"].str.split(":", expand=True)
    df[["start", "end"]] = df["coordinates"].str.split("\\.\\.", expand=True)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df = df.sort_values(by=["scaffold", "strand", "start"]).copy()
    df.drop(columns=["coordinates"], inplace=True)
else:
    df["scaffold"] = df["seqName"]


def get_dev_condition(lib_name):
    if tool_name == "sRNAbench":
        return lib_name.split("_")[1][:-1]
    else:
        return lib_name[:-1]

# Lists to store results
good_candidates = []
filter_out = []

# Iterate through each (scaffold, strand) to build groups
for (scaffold, strand), sub_df in df.groupby(["scaffold", "strand"]):
    sub_df = sub_df.copy()
    sub_df["dev_condition"] = sub_df[dev_condition_col].apply(get_dev_condition)
    current_group = []
    first_start = None

    for _, row in sub_df.iterrows():
        if not current_group:
            first_start = row["start"]
            current_group.append(row)
        elif row["start"] <= first_start + 20:
            current_group.append(row)
        else:
            if len(current_group) > 1:
                dev_counts = pd.Series([x["dev_condition"] for x in current_group]).value_counts()
                if any(dev_counts >= 2):
                    highest_read_seq = max(current_group, key=lambda x: x[read_count_col])
                    same_precursor = len(set(x[precursor_col] for x in current_group)) == 1
                    same_mature_star = len(set((x[mature_col], x[star_col]) for x in current_group)) == 1
                    highest_read_seq["all_same"] = "yes" if same_precursor and same_mature_star else "no"
                    highest_read_seq["overlaps"] = len(current_group)
                    good_candidates.append(highest_read_seq)
                else:
                    filter_out.append(current_group[0])
            current_group = [row]
            first_start = row["start"]

    if len(current_group) > 1:
        dev_counts = pd.Series([x["dev_condition"] for x in current_group]).value_counts()
        if any(dev_counts >= 2):
            highest_read_seq = max(current_group, key=lambda x: x[read_count_col])
            same_precursor = len(set(x[precursor_col] for x in current_group)) == 1
            same_mature_star = len(set((x[mature_col], x[star_col]) for x in current_group)) == 1
            highest_read_seq["all_same"] = "yes" if same_precursor and same_mature_star else "no"
            highest_read_seq["overlaps"] = len(current_group)
            good_candidates.append(highest_read_seq)
        else:
            filter_out.append(current_group[0])

# Convert lists to DataFrames
good_candidates_df = pd.DataFrame(good_candidates)
filter_out_df = pd.DataFrame(filter_out)

# Save output files
good_candidates_df.to_csv(output_path + f"{tool_name}_goodCandidates.csv", index=False)
filter_out_df.to_csv(output_path + f"{tool_name}_filterout.csv", index=False)

print(f"Processing complete for {tool_name}. Check output files.")
