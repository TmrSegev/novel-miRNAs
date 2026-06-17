#!/usr/bin/python
"""Deprecated wrapper — use srnabenchPerLibraryFilter.py."""
from pipeline_wrapper import run_canonical

if __name__ == "__main__":
    run_canonical(
        "srnabenchPerLibraryFilter.py",
        message="hofsteniasRNAbenchFilter.py is deprecated; use srnabenchPerLibraryFilter.py",
    )
