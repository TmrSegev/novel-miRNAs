#!/usr/bin/python
"""Deprecated wrapper — use srnabenchUniteGFF.py."""
from pipeline_wrapper import run_canonical

if __name__ == "__main__":
    run_canonical(
        "srnabenchUniteGFF.py",
        message="nematodesRNAbenchGFF.py is deprecated; use srnabenchUniteGFF.py",
    )
