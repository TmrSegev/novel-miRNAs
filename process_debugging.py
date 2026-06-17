#!/usr/bin/python
"""Deprecated wrapper — use processGoodCandidates.py."""
from pipeline_wrapper import run_canonical

if __name__ == "__main__":
    run_canonical(
        "processGoodCandidates.py",
        message="process_debugging.py is deprecated; use processGoodCandidates.py",
    )
