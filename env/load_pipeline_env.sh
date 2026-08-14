#!/usr/bin/env bash
# Sourceable cluster env for novel-miRNAs manual runs.
#
# Usage (MUST be sourced, not executed):
#   source /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/env/load_pipeline_env.sh Macrosperma
#   source .../env/load_pipeline_env.sh Macrosperma new_genome
#   source .../env/load_pipeline_env.sh Hofstenia_newGenome
#
# After sourcing: SPECIES, TRACK, SPECIES_DIR, LIBRARIES, verify helpers, etc.
# Keep species path tables in sync with docs/pipeline_v4.md + pipeline_config.py.

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "ERROR: source this file, do not execute it:" >&2
  echo "  source ${BASH_SOURCE[0]} <Species|Species_newGenome> [new_genome]" >&2
  exit 1
fi

_nm_usage() {
  cat >&2 <<'EOF'
Usage:
  source env/load_pipeline_env.sh <Species> [new_genome]
  source env/load_pipeline_env.sh <Species_newGenome>

Species: Elegans | Macrosperma | Sulstoni | Hofstenia
EOF
}

_NM_ARG1="${1:-}"
_NM_ARG2="${2:-}"
if [[ -z "$_NM_ARG1" ]]; then
  _nm_usage
  return 1 2>/dev/null || exit 1
fi

# Resolve script dir → REPO (works when sourced from any cwd)
_NM_ENV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO="$(cd "$_NM_ENV_DIR/.." && pwd)"
export BASE="${BASE:-/mnt/new_groups/vaksler_group/Isana_Tzah/Charles_seq}"
export RNACENTRAL="${RNACENTRAL:-/mnt/new_groups/vaksler_group/Isana_Tzah/RNAcentral}"

# Accept Species_newGenome as a single token
if [[ "$_NM_ARG1" == *_newGenome ]]; then
  export SPECIES="${_NM_ARG1%_newGenome}"
  export VARIANT="--variant new_genome"
elif [[ "$_NM_ARG2" == "new_genome" || "$_NM_ARG2" == "--variant" ]]; then
  export SPECIES="$_NM_ARG1"
  export VARIANT="--variant new_genome"
else
  export SPECIES="$_NM_ARG1"
  export VARIANT=""
fi

case "$SPECIES" in
  Elegans|Macrosperma|Sulstoni|Hofstenia) ;;
  *)
    echo "ERROR: unknown species '$SPECIES'" >&2
    _nm_usage
    return 1 2>/dev/null || exit 1
    ;;
esac

if [[ -n "$VARIANT" ]]; then
  export TRACK="${SPECIES}_newGenome"
else
  export TRACK="$SPECIES"
fi

# Libraries from pipeline_config.py (stays in sync with code)
if ! export LIBRARIES="$(
  cd "$REPO" && python3 - <<PY
from pipeline_config import SPECIES_CONFIG
print(",".join(SPECIES_CONFIG["${SPECIES}"]["libraries"]))
PY
)"; then
  echo "ERROR: could not read LIBRARIES from pipeline_config.py" >&2
  return 1 2>/dev/null || exit 1
fi

export SPECIES_DIR="$BASE/$TRACK"
export SCRIPTS_DIR="$SPECIES_DIR/scripts"
export RNA_MI_DIR="$RNACENTRAL/miRNAs/$TRACK"
export BLAST_QUERY_DIR="$RNACENTRAL/queries/$TRACK"
export ZIV_XLSX="$BASE/Ziv_Features/all_remaining_after_ziv_${TRACK}.xlsx"
export INTERSECTIONS_XLSX="$RNA_MI_DIR/intersections_table_${SPECIES}.xlsx"
export SEED="$BASE/mirbase_data/Seeds.txt"

# Species-specific path layout (old vs new_genome)
case "$SPECIES" in
  Elegans)
    export INDEX_BASENAME=elegans
    export HOF_FLAGS=""
    export ZIV_SHEET="(D) Structural Features"
    export MIRGE_FASTA_DIR="$SPECIES_DIR/miRge/"
    if [[ -n "$VARIANT" ]]; then
      export BASH_DIR="$SPECIES_DIR/bash"
      export GENOME_DIR="$SPECIES_DIR/genome"
      export GENOME_FA="$GENOME_DIR/CELEG.caenorhabditis_elegans_PRJNA13758_WBPS19.scaffolds.fna"
      export GENOME_FA_NO_WS="$GENOME_FA"
      export SRNABENCH_INDEX=elegansNewGenomeIndexed
    else
      export BASH_DIR="$SPECIES_DIR/Bash"
      export GENOME_DIR="$SPECIES_DIR/Genome"
      export GENOME_FA="$GENOME_DIR/caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa"
      export GENOME_FA_NO_WS="$GENOME_DIR/new_caenorhabditis_elegans.PRJNA13758.WBPS16.genomic.fa"
      export SRNABENCH_INDEX=elegansGenomeIndexed
    fi
    unset READ_FASTQ_DIR
    ;;
  Macrosperma)
    export INDEX_BASENAME=macrosperma
    export HOF_FLAGS=""
    export ZIV_SHEET="(D) Structural Features"
    export MIRGE_FASTA_DIR="$SPECIES_DIR/miRge/"
    export BASH_DIR="$SPECIES_DIR/bash"
    export GENOME_DIR="$SPECIES_DIR/genome"
    if [[ -n "$VARIANT" ]]; then
      export GENOME_FA="$GENOME_DIR/CMACR.caenorhabditis_macrosperma_JU2083_v2.scaffolds.fna"
      export SRNABENCH_INDEX=macrospermaNewGenomeIndexed
    else
      export GENOME_FA="$GENOME_DIR/CMACR.caenorhabditis_macrosperma_JU2083_v1.scaffolds.fna"
      export SRNABENCH_INDEX=macrospermaGenomeIndexed
    fi
    export GENOME_FA_NO_WS="$GENOME_FA"
    unset READ_FASTQ_DIR
    ;;
  Sulstoni)
    export INDEX_BASENAME=sulstoni
    export HOF_FLAGS=""
    export ZIV_SHEET="(D) Structural Features"
    export MIRGE_FASTA_DIR="$SPECIES_DIR/miRge/"
    export BASH_DIR="$SPECIES_DIR/bash"
    export GENOME_DIR="$SPECIES_DIR/genome"
    if [[ -n "$VARIANT" ]]; then
      export GENOME_FA="$GENOME_DIR/CSULS.caenorhabditis_sulstoni_PRJEB12601_WBPS19.scaffolds.fna"
      export SRNABENCH_INDEX=sulstoniNewGenomeIndexed
    else
      export GENOME_FA="$GENOME_DIR/CSULS.caenorhabditis_sulstoni_JU2788_v1.scaffolds.fna"
      export SRNABENCH_INDEX=sulstoniGenomeIndexed
    fi
    export GENOME_FA_NO_WS="$GENOME_FA"
    unset READ_FASTQ_DIR
    ;;
  Hofstenia)
    export INDEX_BASENAME=hofstenia
    export HOF_FLAGS="--base-path $BASE"
    export ZIV_SHEET="(A) Unfiltered"
    export MIRGE_FASTA_DIR="$SPECIES_DIR/miRge_after_Ziv/"
    export BASH_DIR="$SPECIES_DIR/bash"
    # Reads always from the original Hofstenia Fastq tree
    export READ_FASTQ_DIR="$BASE/Hofstenia/Fastq/Hmia_annotation/filtered"
    if [[ -n "$VARIANT" ]]; then
      export GENOME_DIR="$SPECIES_DIR/sRNA_PBonly"
      export GENOME_FA="$GENOME_DIR/hofPB_v6.FINAL.fa"
      export SRNABENCH_INDEX=hofsteniaNewGenomeIndexed
    else
      export GENOME_DIR="$SPECIES_DIR/Genome/refs/Hmia_ref"
      export GENOME_FA="$GENOME_DIR/Hmia.030120.fasta"
      export SRNABENCH_INDEX=hofsteniaGenomeIndexed
    fi
    export GENOME_FA_NO_WS="$GENOME_FA"
    ;;
esac

# Always overwrite so a previous `nm` cannot leak SEQOBJ_NAME across tracks.
# Hofstenia_newGenome sRNAbench wrappers use species=hofPB_v6 (zip from makeSeqObj).
if [[ "$TRACK" == "Hofstenia_newGenome" ]]; then
  export SEQOBJ_NAME=hofPB_v6
else
  export SEQOBJ_NAME="$SRNABENCH_INDEX"
fi
export SEQOBJ_ZIP="$BASE/sRNAtoolboxDB/seqOBJ/${SEQOBJ_NAME}.zip"

export STAR_SAMS="$(for lib in ${LIBRARIES//,/ }; do echo ../STAR/align_to_genome/$lib/${SPECIES}_Aligned.out.sam; done)"

# Verify helpers (same as docs/pipeline_v4.md)
FAIL=0
ok()   { echo "OK: $*"; }
fail() { echo "FAIL: $*"; FAIL=1; }
# Pipeline outputs: must exist, be non-empty, and have mtime within the last 7 days.
need_file() {
  local f="$1"
  if [[ ! -s "$f" ]]; then
    fail "missing/empty: $f"
  elif [[ -z "$(find "$f" -mtime -7 2>/dev/null)" ]]; then
    fail "stale (mtime >7d): $f"
  else
    ok "file $f ($(wc -c <"$f") bytes, mtime≤7d)"
  fi
}
# Reference inputs (e.g. genome FASTA): existence only — not produced by this phase.
need_input() {
  local f="$1"
  if [[ -s "$f" ]]; then ok "input $f ($(wc -c <"$f") bytes)"
  else fail "missing/empty input: $f"; fi
}
need_dir() {
  local d="$1"
  if [[ -d "$d" ]]; then ok "dir $d"
  else fail "missing dir: $d"; fi
}
count_lines() { wc -l <"$1" | tr -d ' '; }

# Convenience: snapshot current track
nm_snapshot() {
  local snap="${BASE}/snapshots/${TRACK}_$(date +%Y%m%d_%H%M%S)"
  mkdir -p "$snap"
  [[ -d "$SCRIPTS_DIR" ]] && cp -a "$SCRIPTS_DIR" "$snap/scripts"
  [[ -d "$SPECIES_DIR/unique_candidates" ]] && cp -a "$SPECIES_DIR/unique_candidates" "$snap/unique_candidates"
  [[ -d "$SPECIES_DIR/counts_sep" ]] && cp -a "$SPECIES_DIR/counts_sep" "$snap/counts_sep"
  [[ -d "$RNA_MI_DIR" ]] && cp -a "$RNA_MI_DIR" "$snap/miRNAs_${TRACK}"
  [[ -f "$ZIV_XLSX" ]] && cp -a "$ZIV_XLSX" "$snap/"
  [[ -d "$BLAST_QUERY_DIR" ]] && cp -a "$BLAST_QUERY_DIR" "$snap/blast_queries"
  echo "Snapshot at $snap"
}

echo "Loaded pipeline env:"
echo "  SPECIES=$SPECIES  TRACK=$TRACK  VARIANT=${VARIANT:-<empty>}"
echo "  SPECIES_DIR=$SPECIES_DIR"
echo "  LIBRARIES=$LIBRARIES"
echo "  GENOME_FA=$GENOME_FA"
echo "  SRNABENCH_INDEX=$SRNABENCH_INDEX  SEQOBJ_NAME=$SEQOBJ_NAME"
echo "  SEQOBJ_ZIP=$SEQOBJ_ZIP"
echo "Helpers: need_file (≤7d) / need_input / need_dir / nm_snapshot"
unset _NM_ARG1 _NM_ARG2 _NM_ENV_DIR
