# novel-miRNAs — agent notes

## How this repo is operated

- **Edit locally** on the laptop, **push to git**, then **pull on the cluster**.
- Cluster path for scripts: `/mnt/new_groups/vaksler_group/Isana_Tzah/novel-miRNAs`
- Data and job submission live on the cluster; coding agents do not run the pipeline there.
- Every command in `docs/pipeline_v3.md` must be **copy-paste ready** for a human SSH session (exported vars + absolute `$REPO/...` invocations). Do not assume an agent will execute those phases.
