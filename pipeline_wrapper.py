"""Helpers for backward-compatible script wrappers."""

import os
import subprocess
import sys
import warnings


def translate_argv(argv):
    """Translate legacy flags (e.g. --new-genome) to canonical form."""
    out = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg == "--new-genome":
            nxt = argv[i + 1] if i + 1 < len(argv) else None
            if nxt in (None, "True", "true"):
                out.extend(["--variant", "new_genome"])
                i += 2 if nxt is not None else 1
                continue
            if nxt in ("False", "false"):
                i += 2
                continue
            out.extend(["--variant", "new_genome"])
            i += 1
            continue
        out.append(arg)
        i += 1
    return out


def run_canonical(canonical_name, argv=None, message=None):
    """Re-exec the canonical script with translated argv."""
    argv = translate_argv(argv if argv is not None else sys.argv[1:])
    here = os.path.dirname(os.path.abspath(__file__))
    target = os.path.join(here, canonical_name)
    if message:
        warnings.warn(message, DeprecationWarning, stacklevel=2)
    cmd = [sys.executable, target] + argv
    sys.exit(subprocess.call(cmd))
