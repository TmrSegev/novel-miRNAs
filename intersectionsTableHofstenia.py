#!/usr/bin/python
"""Deprecated wrapper — use intersectionsTable.py."""
from pipeline_wrapper import run_canonical

if __name__ == "__main__":
    run_canonical(
        "intersectionsTable.py",
        message="intersectionsTableHofstenia.py is deprecated; use intersectionsTable.py",
    )
