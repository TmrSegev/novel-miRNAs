#!/usr/bin/python
"""Deprecated wrapper — use mirdeepPerLibraryFilter.py."""
from pipeline_wrapper import run_canonical

if __name__ == "__main__":
    run_canonical(
        "mirdeepPerLibraryFilter.py",
        message="hofsteniaMirdeepFilter.py is deprecated; use mirdeepPerLibraryFilter.py",
    )
