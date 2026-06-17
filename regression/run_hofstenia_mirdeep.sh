#!/usr/bin/env bash
# Sandbox regression: Hofstenia miRDeep unite + good_candidates (pass 1 + 2).
# Reads production inputs via symlinks; writes only under regression/sandbox/.

set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
PROD="/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq"
SANDBOX="$REPO/regression/sandbox/Charles_seq"
PROD_SCRIPTS="$PROD/Hofstenia/scripts"

rm -rf "$SANDBOX/Hofstenia/scripts" "$SANDBOX/Hofstenia/good_candidates"
mkdir -p "$SANDBOX/Hofstenia/scripts" "$SANDBOX/Hofstenia/good_candidates"
ln -sfn "$PROD/Hofstenia/mirdeep_out" "$SANDBOX/Hofstenia/mirdeep_out"
ln -sfn "$PROD/mirbase_data" "$SANDBOX/mirbase_data"

cd "$REPO"

echo "=== Pass 1: debugging CSV ==="
python3 mirdeepUniteGFF.py -o Hofstenia_mirdeep.gff3 --create-fasta Hofstenia_mirdeep.fasta \
  -s Hofstenia --base-path "$SANDBOX" --goodcandidates False

echo "=== processGoodCandidates ==="
python3 processGoodCandidates.py --tool miRDeep -s Hofstenia --base-path "$SANDBOX"

echo "=== Pass 2: final GFF ==="
python3 mirdeepUniteGFF.py -o Hofstenia_mirdeep.gff3 --create-fasta Hofstenia_mirdeep.fasta \
  -s Hofstenia --base-path "$SANDBOX" --goodcandidates True

echo "=== Diff vs production ==="
for f in debugging_Hofstenia_miRDeep_1.csv debugging_Hofstenia_miRDeep_2.csv Hofstenia_mirdeep.gff3; do
  echo "--- $f ---"
  diff -q "$PROD_SCRIPTS/$f" "$SANDBOX/Hofstenia/scripts/$f" && echo IDENTICAL || echo DIFFERS
done
