#!/usr/bin/env bash
# Optional one-line cluster ~/.bashrc hook for MobaXterm / SSH.
#
# Add to ~/.bashrc on the cluster (once):
#   source /mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs/env/moba_aliases.sh
#
# Then in any new session:
#   nm Macrosperma
#   nm Macrosperma new_genome
#   nm Hofstenia_newGenome
#   nm-list

_NM_REPO="${NM_REPO:-/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs}"
_NM_LOADER="$_NM_REPO/env/load_pipeline_env.sh"

nm() {
  if [[ ! -f "$_NM_LOADER" ]]; then
    echo "ERROR: missing $_NM_LOADER (git pull in novel-miRNAs?)" >&2
    return 1
  fi
  # shellcheck disable=SC1090
  source "$_NM_LOADER" "$@"
}

nm-list() {
  cat <<'EOF'
Tracks:
  nm Elegans
  nm Elegans new_genome          # or: nm Elegans_newGenome
  nm Macrosperma
  nm Macrosperma new_genome
  nm Sulstoni
  nm Sulstoni new_genome
  nm Hofstenia
  nm Hofstenia new_genome

After load:
  echo $SPECIES_DIR
  nm_snapshot          # backup current track outputs
EOF
}

echo "novel-miRNAs aliases ready (nm, nm-list). Repo: $_NM_REPO"
