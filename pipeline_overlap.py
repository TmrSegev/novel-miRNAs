"""Coordinate overlap deduplication for united candidate tables."""

import pandas as pd
import numpy as np


def deduplicate_coordinate_overlaps(table, scaffold_col, min_fraction=0.6):
    """
    Collapse rows with >= min_fraction coordinate overlap on the same scaffold.
    Keeps one representative per cluster; records overlap count in 'overlaps'.
    Rows with no overlapping partner are moved to removed_no_overlaps.
    """
    table = table.sort_values([scaffold_col, "start", "end"]).copy()
    table["overlaps"] = np.zeros(len(table), dtype=int)
    no_overlaps = pd.DataFrame(columns=table.columns)

    for index, row in table.iterrows():
        if index not in table.index:
            continue
        span = row["end"] - row["start"]
        if span <= 0:
            continue
        table["distance"] = (row["end"] - table["start"]) / span
        overlaps = table[(table["distance"] >= min_fraction) & (table["distance"] <= 1)].tail(-1)
        overlaps = overlaps[overlaps[scaffold_col] == row[scaffold_col]]
        table.loc[index, "overlaps"] = len(overlaps)
        if len(overlaps) == 0:
            no_overlaps = pd.concat([no_overlaps, pd.DataFrame([row])], ignore_index=True)
            table = table.drop(index)
        else:
            table = table.drop(overlaps.index)

    if "distance" in table.columns:
        table = table.drop(columns=["distance"])
    return table, no_overlaps
