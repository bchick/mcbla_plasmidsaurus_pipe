# TODO

Goal: this pipeline reproduces the Plasmidsaurus portal outputs (BAMs and
count tables) from raw FASTQ, and others in the lab can run it easily.
Verified recipe and results: `validation/RESULTS_2026-09-24.md`.

## 1. Port the verified recipe into the main pipeline

- [x] `scripts/01_filter_reads.sh`: fastp `--cut_tail` mean Q20 (was Q15)
- [x] `scripts/02_align.sh`: STAR defaults + `--outFilterIntronMotifs RemoveNoncanonical`;
      dropped BySJout, mismatch ratio 0.04, SJ overhang 8, intron min/max,
      multimap 20; added `NM`; unsorted output, sorted by samtools in step 3;
      reads written in FASTQ order (`--outSAMorder PairedKeepInputOrder`) so
      dedup no longer changes between runs
- [x] `scripts/04_dedup.sh`: UMI separator `_`; UMICollapse flags fixed
      (`--algo dir -k 1`; the old `--algo directional --edit-distance` crashed it);
      Java heap 64g (the wrapper default of 4g ran out of memory)
- [x] `scripts/07_quantify.sh`: `-t gene`, MAPQ 0, strand default 1 (also
      `run_pipeline.sh` fallback)
- [x] Configs: portal recipe values, `paired_end: false`; mouse config points
      at Ensembl 114 + ERCC in `data/reference/` (relative paths resolved
      against the repo root); BED12 for RSeQC made with UCSC gtfToGenePred
- [ ] Human reference: build Ensembl GRCh38 + ERCC and update `config.human.yaml`
      (still GENCODE paths)
- [x] `scripts/05_mapping_qc.sh`: gene body coverage on housekeeping genes
      (`scripts/make_housekeeping_bed.py`; ~10 min/sample instead of hours);
      RSeQC `log.txt` no longer written to the launch directory
- [x] Acceptance test (`validation/compare_to_portal.py`): PASS on all 4
      LJQQSK samples (2026-09-26, `validation/ACCEPTANCE_2026-09-26.txt`).
      STAR exact vs portal; ~0.1% of genes differ by at most 3 reads
      (read order); 32- and 16-thread runs identical.
  - [x] Confirm two runs (32 and 16 threads) give identical results
  - [x] Revise pass criteria: exact vs STAR stats and vs a repeat run; small
        tolerance vs the portal, whose read order cannot be reproduced
  - [x] Update `validation/RESULTS_2026-09-24.md` with the read-order finding
- [ ] Run steps 8-10 (edgeR, GSEA) on the final counts; untested since the
      GSEA fixes
- [ ] Check QC steps 5-6 finished on the full run (housekeeping gene body
      coverage, MultiQC)

## 2. Make it easy for others to use

- [ ] One-command setup via pixi (`pixi install`, then `pixi run pipeline ...`)
- [ ] Scripted reference build (Ensembl GRCm39 r114 + ERCC92 + BED12; human equivalent;
      add UCSC gtfToGenePred/genePredToBed to pixi)
- [ ] Automatic single-end vs paired-end detection
- [ ] Package the portal-concordance check as a test anyone can run
- [ ] README: quickstart, inputs/outputs, how the output relates to the portal.
      Also fix stale recipe text in README and CLAUDE.md (e.g. "exons and
      3' UTR" counting; BCL Convert/fqtk listed as step 1)
- [ ] Sample IDs keep the portal's `_r1`/`_r2` suffix (e.g.
      `LJQQSK_1_BR4_1_US_r1`); check `samples.tsv` naming matches

## 3. Housekeeping

- [x] Commit today's work on `fix/portal-concordance`
- [x] Track only scripts and write-ups under `validation/`; run outputs
      (`LJQQSK/`, `work/`) are ignored
- [x] Add `__pycache__/` to `.gitignore`
- [ ] Root filesystem (`/`, holding `/tmp`) at 88%: keep big intermediates on `/data`
- [ ] Delete `validation/work/` once the acceptance test is final (~20 GB)

## Low priority

- [x] Dedup residual (3-28 multi-mapped reads per sample): read order at tied
      positions (STAR threading) decides which duplicate UMICollapse keeps.
      Our runs are now reproducible; the portal's order cannot be recovered.
