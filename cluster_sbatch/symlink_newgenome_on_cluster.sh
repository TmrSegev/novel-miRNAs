#!/bin/bash
# Run ON THE CLUSTER after: cd "$REPO" && git pull
# Makes Charles_seq / RNAcentral new-genome sbatch paths point at cluster_sbatch/.
set -euo pipefail

REPO="${REPO:-/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs}"
BASE="${BASE:-/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq}"
RNAC="${RNACENTRAL:-/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral}"
TS=$(date +%Y%m%d_%H%M%S)
ARCH="${BASE}/_archive/newgenome_sbatch_before_symlink_${TS}"
mkdir -p "$ARCH"

link_file() {
  local src="$1" dest="$2"
  mkdir -p "$(dirname "$dest")"
  if [[ -e "$dest" && ! -L "$dest" ]]; then
    local rel="${dest#$BASE/}"
    rel="${rel#$RNAC/}"
    mkdir -p "$ARCH/$(dirname "$rel")"
    mv "$dest" "$ARCH/$rel"
  fi
  ln -sfn "$src" "$dest"
  echo "  link $dest"
}

echo "Archive: $ARCH"

# --- nematode bash: whole-directory symlink (thin) ---
for sp in Elegans Macrosperma Sulstoni; do
  dest="$BASE/${sp}_newGenome/bash"
  src="$REPO/cluster_sbatch/${sp}_newGenome"
  echo "=== $sp bash ==="
  if [[ -L "$dest" ]]; then
    ln -sfn "$src" "$dest"
    echo "  refreshed symlink $dest"
  else
    mkdir -p "$(dirname "$dest")"
    if [[ -e "$dest" ]]; then
      mv "$dest" "$ARCH/${sp}_newGenome_bash"
    fi
    ln -sfn "$src" "$dest"
    echo "  $dest -> $src"
  fi
done

# --- Hofstenia bash: per-file only (mapper run dirs live here) ---
echo "=== Hofstenia bash (per-file) ==="
for f in "$REPO/cluster_sbatch/Hofstenia_newGenome/"*.sbatch; do
  bn=$(basename "$f")
  link_file "$f" "$BASE/Hofstenia_newGenome/bash/$bn"
done

# --- filter sbatch in scripts/ (scripts/ also holds GFF/CSV) ---
for sp in Elegans Macrosperma Sulstoni Hofstenia; do
  echo "=== $sp scripts filters ==="
  for f in "$REPO/cluster_sbatch/scripts/${sp}_newGenome/"*.sbatch; do
    [[ -f "$f" ]] || continue
    bn=$(basename "$f")
    link_file "$f" "$BASE/${sp}_newGenome/scripts/$bn"
  done
done

# --- per-library mirdeep_test ---
for sp in Elegans Macrosperma Sulstoni Hofstenia; do
  echo "=== $sp mirdeep_test ==="
  for d in "$REPO/cluster_sbatch/mirdeep_test/${sp}_newGenome/"*; do
    [[ -d "$d" ]] || continue
    lib=$(basename "$d")
    src="$d/mirdeep_test.sbatch"
    [[ -f "$src" ]] || continue
    link_file "$src" "$BASE/${sp}_newGenome/mirdeep_out/$lib/mirdeep_test.sbatch"
  done
done

# --- RNAcentral intersections ---
for sp in Elegans Macrosperma Sulstoni Hofstenia; do
  echo "=== RNAcentral $sp newGenome intersections ==="
  srcdir="$REPO/cluster_sbatch/RNAcentral/miRNAs/${sp}_newGenome"
  destdir="$RNAC/miRNAs/${sp}_newGenome"
  mkdir -p "$destdir"
  for f in "$srcdir/"*.sbatch; do
    [[ -f "$f" ]] || continue
    link_file "$f" "$destdir/$(basename "$f")"
  done
done

# --- BLAST new-genome query wrappers (RNAcentral/bash is a real dir) ---
echo "=== BLAST newGenome query wrappers ==="
for f in "$REPO/cluster_sbatch/RNAcentral/bash/"*_newgenome_queries.sbatch; do
  [[ -f "$f" ]] || continue
  link_file "$f" "$RNAC/bash/$(basename "$f")"
done

echo
echo "Done. Spot-check:"
ls -ld "$BASE/Elegans_newGenome/bash" "$BASE/Hofstenia_newGenome/bash"
ls -l "$BASE/Elegans_newGenome/bash/featurecounts_mirdeep_sep.sbatch"
ls -l "$BASE/Hofstenia_newGenome/bash/sRNAbench_EC1.sbatch"
ls -l "$RNAC/bash/blast_elegans_newgenome_queries.sbatch"
