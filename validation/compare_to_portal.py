#!/usr/bin/env python3
"""
compare_to_portal.py - Acceptance test: run_pipeline.sh output vs the portal

Description:
    Compares the output of run_pipeline.sh (steps 1-7) on a Plasmidsaurus run
    against the processed files the Plasmidsaurus portal delivered for the same
    run: STAR mapping statistics, deduplicated BAMs, and the gene count matrix.
    Optionally also compares against a previous validation/reproduce_portal.sh
    run, which used the verified recipe directly and should agree exactly.

    This is the acceptance test for porting the verified recipe into the main
    pipeline (see validation/RESULTS_2026-09-24.md).

Usage:
    python3 validation/compare_to_portal.py \
        --pipeline-out validation/work/acceptance/out \
        --portal data/plasmidsaurus/LJQQSK \
        --run LJQQSK \
        [--recipe-out validation/LJQQSK]

Inputs:
    --pipeline-out  run_pipeline.sh output directory (02_aligned/, 04_dedup/,
                    07_counts/gene_counts.txt)
    --portal        Portal download: results/<RUN>-expression-matrix.tsv,
                    results/<RUN>-mapping-stats-reads.csv, <RUN>_bam/*.bam
    --run           Plasmidsaurus run ID (file name prefix)
    --recipe-out    Optional reproduce_portal.sh output (counts.txt, <sid>/dedup.bam)

Outputs:
    Report on stdout. Exit status 0 if every check passes, 1 otherwise.

Pass criteria (what the verified recipe achieved on LJQQSK):
    - STAR uniquely mapped and multi-mapped read counts equal the portal's
    - Every gene within 1 read of the portal's count; log1p Pearson r > 0.9999
    - If --recipe-out is given: count matrix identical to the recipe run

Dependencies:
    Python >= 3.8 with pandas, numpy; samtools on PATH (use `pixi run`)

Author: Brent Chick
Date: 2026-09-26
Version: 1.0.0
"""

import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd


# ==============================================================================
# HELPERS
# ==============================================================================

def run_prefix(sample_id):
    """Short portal sample ID, e.g. 'LJQQSK_1' from 'LJQQSK_1_BR4_1_US_r1'."""
    return "_".join(sample_id.split("_")[:2])


def primary_nh_counts(bam):
    """Primary mapped reads in a BAM, split into unique (NH==1) and multi (NH>1)."""
    cmd = f"samtools view -F 260 '{bam}' | grep -o 'NH:i:[0-9]*' | sort | uniq -c"
    out = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True).stdout
    counts = {int(nh.split(":")[-1]): int(n)
              for n, nh in (line.split() for line in out.strip().splitlines())}
    return counts.get(1, 0), sum(v for k, v in counts.items() if k > 1)


def star_stat(log_text, key):
    """Integer value of a field in STAR's Log.final.out."""
    return int(re.search(rf"{key} \|\s+(\d+)", log_text).group(1))


def read_featurecounts(path):
    """featureCounts matrix as genes x samples, columns named by sample ID.

    Keeps only the count columns (those that were BAM paths), so it works
    whether or not --extraAttributes added annotation columns.
    """
    df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    df = df[[c for c in df.columns if c.endswith(".bam")]]
    df.columns = [os.path.basename(c) if not c.endswith("/dedup.bam")
                  else c.split("/")[-2] for c in df.columns]
    df.columns = [re.sub(r"(\.dedup)?\.bam$", "", c) for c in df.columns]
    return df


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--portal", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--recipe-out")
    args = ap.parse_args()

    out, portal, run = args.pipeline_out, args.portal, args.run
    stats = pd.read_csv(f"{portal}/results/{run}-mapping-stats-reads.csv", index_col=0)
    ref = pd.read_csv(f"{portal}/results/{run}-expression-matrix.tsv", sep="\t", index_col=0)
    ours = read_featurecounts(f"{out}/07_counts/gene_counts.txt")
    # The portal matrix has every annotated gene; ours should too, but
    # reindex so a missing gene shows up as a count difference, not a crash
    ours = ours.reindex(ref.index).fillna(0)

    recipe = None
    if args.recipe_out:
        recipe = read_featurecounts(f"{args.recipe_out}/counts.txt").reindex(ref.index).fillna(0)

    failures = []

    def check(ok, label):
        if not ok:
            failures.append(label)
        return "PASS" if ok else "FAIL"

    def pct(a, b):
        return f"{(a - b) / b * 100:+.3f}%" if b else "n/a"

    print(f"=== {run}: run_pipeline.sh output vs portal ===")
    print(f"pipeline output: {out}")
    for sid in sorted(ours.columns):
        short = run_prefix(sid)
        # Portal mapping-stats rows are labelled without the run prefix
        label = sid[len(short) + 1:]
        s = stats.loc[label]
        log_text = open(f"{out}/02_aligned/{sid}_Log.final.out").read()
        uniq = star_stat(log_text, "Uniquely mapped reads number")
        multi = star_stat(log_text, "Number of reads mapped to multiple loci")

        pbam = glob.glob(f"{portal}/{run}_bam/{sid}_dedup-mapped-reads.bam")[0]
        pu, pm = primary_nh_counts(pbam)
        ou, om = primary_nh_counts(f"{out}/04_dedup/{sid}.dedup.bam")

        c, r = ours[sid], ref[f"{short}_count"]
        diff = (c - r).abs()
        within1 = np.mean(diff <= 1) * 100
        pearson = np.corrcoef(np.log1p(c), np.log1p(r))[0, 1]

        print(f"\n{sid}")
        print(f"  STAR unique    ours {uniq:>10}  portal {s['Uniquely Mapped']:>10}  "
              f"{pct(uniq, s['Uniquely Mapped'])}  {check(uniq == s['Uniquely Mapped'], f'{sid} STAR unique')}")
        print(f"  STAR multi     ours {multi:>10}  portal {s['Multi-mapped']:>10}  "
              f"{pct(multi, s['Multi-mapped'])}  {check(multi == s['Multi-mapped'], f'{sid} STAR multi')}")
        print(f"  dedup unique   ours {ou:>10}  portal {pu:>10}  {ou - pu:+d} reads")
        print(f"  dedup multi    ours {om:>10}  portal {pm:>10}  {om - pm:+d} reads")
        print(f"  counts total   ours {c.sum():>10.0f}  portal {r.sum():>10.0f}  {c.sum() - r.sum():+.1f}")
        print(f"  genes exact {np.mean(diff < 1e-6) * 100:.2f}%  within 1 read {within1:.2f}% "
              f"{check(within1 == 100, f'{sid} within 1 read')}  "
              f"Pearson r (log1p) {pearson:.6f} {check(pearson > 0.9999, f'{sid} Pearson')}")

        if recipe is not None:
            rmax = (c - recipe[short]).abs().max()
            print(f"  vs recipe run: max |count diff| {rmax:.4g}  "
                  f"{check(rmax < 1e-6, f'{sid} identical to recipe run')}")

    print()
    if failures:
        print(f"RESULT: FAIL ({len(failures)} checks)")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("RESULT: PASS (all checks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
