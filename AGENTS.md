# novel-miRNAs — agent notes

## How this repo is operated

- **Edit locally** on the laptop, **push to git**, then **pull on the cluster**.
- Cluster path for scripts: `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs`
- Data and job submission live on the cluster; coding agents do not run the pipeline there.
- Every command in `docs/pipeline_v3.md` must be **copy-paste ready** for a human SSH session (exported vars + absolute `$REPO/...` invocations). Do not assume an agent will execute those phases.

## Pipeline run / overwrite notes

- Use `$TRACK` (`{Species}` vs `{Species}_newGenome`) so RNA-mi / Ziv / BLAST query paths stay separated. See overwrite + Verify sections in `docs/pipeline_v3.md`.
- **Hofstenia_newGenome** already has prior outputs; nematode `_newGenome` tracks have not been run yet.
- Stop after Phase 11 (Ziv) until Phase 12 (Oscar / 5p-het) is implemented; do not run 13–14 early.
- Snapshot existing tracks before re-running in-place phases (esp. old genomes and Hofstenia_newGenome).
