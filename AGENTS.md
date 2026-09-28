# AGENTS.md — running mcbla-plasmidsaurus-pipe for a lab member

Instructions for AI agents (and people) asked to run this pipeline on the
McBla lab server. It processes **Plasmidsaurus RNA-seq** (single-end, UMI)
FASTQs, not plasmid sequencing. The goal is to reproduce the Plasmidsaurus
portal's BAMs and counts, then run edgeR DE and GSEA.

## 1. Before running anything, ask the user

Ask these questions and wait for the answers. **Do not guess.**

1. **Which genome?** Run `pixi run init --list` to show the options:
   - `mm39` (mouse): Ensembl 114 GRCm39 + ERCC92, the same reference as the
     Plasmidsaurus portal.
   - `hg38` (human): GENCODE v44. It has **no ERCC** sequences, so ERCC
     spike-in reads go unmapped. Tell the user this if they added ERCC.

   Older commands use `-g mm10`, which was a mislabel for this same GRCm39
   reference.
2. **Where are the FASTQs** (the folder Plasmidsaurus delivered), and where
   should the project go (under their own `/data/<user>/`, not inside this
   repo or `/data/resource`)?
3. **Sample conditions and replicates**, and which comparisons they want.
   edgeR contrasts are written `"treatment-control"`.

RNA-seq needs no blacklist, IgG control or spike-in choice. If the user asks
for one, explain that it doesn't apply here.

## 2. Write the project config with `init`

```bash
cd /path/to/mcbla-plasmidsaurus-pipe
pixi install
pixi run init --dir /data/<user>/<project> --genome mm39 --fastq-dir /path/to/fastq
```

`init` takes the reference paths (FASTA, GTF, STAR index, RSeQC BEDs) from
`/data/resource/manifest.yaml` and checks they exist. It writes
`<project>/config.yaml`, a copy of the shipped species config pointing at
those paths, and copies an example `samples.tsv`. It exits with status 2 if an
answer is missing; that means go back and ask the user.

**Never build a STAR index or download a genome yourself**; an index takes
about an hour and 26 GB. If the user needs a genome that isn't listed, stop and
tell them. New references go into `/data/resource` through its maintainer.

## 3. Fill in the sample sheet and contrasts

- `<project>/samples.tsv` has the columns `sample_id condition replicate`.
  `sample_id` must match the FASTQ names without `.fastq.gz`, for example
  `LJQQSK_1_BR4_1_US_r1`.
- Put the user's contrasts under `edger.contrasts` in `<project>/config.yaml`.

Show both to the user before running.

## 4. Dry run, then run

```bash
pixi run bash run_pipeline.sh -i /path/to/fastq -o /data/<user>/<project>/results \
    -c /data/<user>/<project>/config.yaml -m /data/<user>/<project>/samples.tsv --dry-run
pixi run bash run_pipeline.sh -i /path/to/fastq -o /data/<user>/<project>/results \
    -c /data/<user>/<project>/config.yaml -m /data/<user>/<project>/samples.tsv -t 16
```

- Always run through `pixi run` so the pinned tool versions are used. Newer
  versions change the results.
- The server is shared and has no Slurm; keep `-t` at 16–32. A run takes
  hours, so start it under `tmux`/`nohup` and tell the user where the logs are
  (`<project>/results/logs/`).
- To resume, use `-s <step>` (steps 1–10; see `./run_pipeline.sh -h`).
- Don't change the fastp, STAR, UMICollapse or featureCounts settings. They
  reproduce the portal exactly (`validation/RESULTS_2026-09-24.md`).
