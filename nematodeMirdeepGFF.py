#!/usr/bin/python
"""Deprecated wrapper — use mirdeepUniteGFF.py."""
from pipeline_wrapper import run_canonical

if __name__ == "__main__":
    run_canonical(
        "mirdeepUniteGFF.py",
        message="nematodeMirdeepGFF.py is deprecated; use mirdeepUniteGFF.py",
    )
