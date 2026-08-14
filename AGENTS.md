# novel-miRNAs — agent notes

## How this repo is operated

- **Edit locally** on the laptop, **push to git**, then **pull on the cluster**.
- **Auto-commit + push:** After each meaningful change in this repo, commit and `git push` to `origin` without asking. Do not wait for an explicit “commit/push” request. Still never force-push, amend others’ commits, or rewrite history unless explicitly asked.
- Cluster path for scripts: `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs`
- Data and job submission live on the cluster; coding agents do not run the pipeline there.
- Every command in `docs/pipeline_v4.md` must be **copy-paste ready** for a human SSH session (exported vars + absolute `$REPO/...` invocations). Do not assume an agent will execute those phases.

## Job wrappers (`cluster_sbatch/`)

- Slurm wrappers live in-repo under `cluster_sbatch/` (source of truth):
  - `{Species}/` — bash/`Bash` wrappers (symlinked from `$BASH_DIR` on the cluster)
  - `{Species}_newGenome/` — featureCounts wrappers for the alternate-assembly track (copy into `$BASH_DIR` after `git pull`; nematode/Hofstenia `_newGenome/bash` is not yet a full-dir symlink)
  - `scripts/{Species}/` — filter `*.sbatch` (symlinked from `$SPECIES_DIR/scripts/`)
  - `mirdeep_test/{Species}/{lib}/` — per-library `mirdeep_test.sbatch`
  - `RNAcentral/` — BLAST + `intersections.sbatch`
- On the cluster, old-genome Charles_seq / RNAcentral paths are **symlinks** into this tree; edit in git, then `git pull` — do not hand-edit the symlink targets as separate copies.
- **Hofstenia is the canonical shared pattern** (per-library discovery, featureCounts naming like `miRNA_miRdeep_*`). Nematodes follow it unless species-specific (path casing, nested sRNAbench out, BLAST, miRBase).

## Pipeline run / overwrite notes

- Use `$TRACK` (`{Species}` vs `{Species}_newGenome`) so RNA-mi / Ziv / BLAST query paths stay separated. See overwrite + Verify sections in `docs/pipeline_v4.md`.
- Cluster env loader (source, do not execute): `env/load_pipeline_env.sh` — or `env/moba_aliases.sh` → `nm Macrosperma`. Libraries come from `pipeline_config.py`.
- **Hofstenia_newGenome** already has prior outputs; nematode `_newGenome` tracks have not been run yet.
- Stop after Phase 11 (Ziv) until Phase 12 (Oscar / 5p-het) is implemented; do not run 13–14 early.
- Snapshot existing tracks before re-running in-place phases (esp. old genomes and Hofstenia_newGenome).
