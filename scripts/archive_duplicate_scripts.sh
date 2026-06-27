#!/usr/bin/env bash
# Rename duplicate pipeline scripts outside novel-miRNAs/ to *_archived.py.
# Idempotent: skips files already archived or missing.
set -euo pipefail

ROOT="/mnt/new_groups/vaksler_group/Isana_Tzah"
REPO="$ROOT/novel-miRNAs"
MANIFEST="$REPO/docs/archived_scripts_manifest.txt"
DRY_RUN="${1:-}"

python3 - "$ROOT" "$REPO" "$MANIFEST" "$DRY_RUN" << 'PY'
import sys
import unicodedata
from datetime import date
from pathlib import Path

root = Path(sys.argv[1])
repo = Path(sys.argv[2])
manifest = Path(sys.argv[3])
dry_run = sys.argv[4] == "--dry-run"

def normalize_basename(name):
    return "".join(c for c in name if unicodedata.category(c) not in ("Cf",))

repo_basenames = {}
for p in repo.glob("*.py"):
    repo_basenames[p.name.lower()] = p.name
    repo_basenames[normalize_basename(p.name).lower()] = p.name

duplicates = []
for top in ("Charles_seq", "RNAcentral"):
    for p in (root / top).rglob("*.py"):
        nb = normalize_basename(p.name).lower()
        if nb in repo_basenames or p.name.lower() in repo_basenames:
            if "_archived" not in p.stem:
                duplicates.append(p)

def archive_name(path: Path) -> Path:
    stem = normalize_basename(path.stem)
    return path.with_name(f"{stem}_archived.py")

lines = [f"# Archived duplicate scripts — {date.today().isoformat()}", ""]
count = 0
for src in sorted(duplicates):
    dst = archive_name(src)
    if dst.exists():
        print(f"SKIP (already archived): {src}")
        lines.append(f"{src} -> {dst}")
        continue
    if not src.exists():
        print(f"SKIP (missing): {src}")
        continue
    print(f"{'DRY-RUN' if dry_run else 'RENAME'}: {src} -> {dst}")
    if not dry_run:
        src.rename(dst)
    lines.append(f"{src} -> {dst}")
    count += 1

if not dry_run:
    manifest.write_text("\n".join(lines) + "\n")
    print(f"\nRenamed {count} files; manifest: {manifest}")
else:
    print(f"\nWould rename {count} files")
PY
