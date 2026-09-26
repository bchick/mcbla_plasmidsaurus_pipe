#!/usr/bin/env python3
"""
compare_to_portal.py - Acceptance test: run_pipeline.sh output vs the portal

Description:
    Compares the output of run_pipeline.sh (steps 1-7) on a Plasmidsaurus run
    against the processed files the Plasmidsaurus portal delivered for the same
    run: STAR mapping statistics, deduplicated BAMs, and the gene count matrix.

    Exact agreement with the portal is not possible after alignment. For reads
    that tie at the same position, UMICollapse keeps the first duplicate in file
    order, and the portal's order came from multithreaded STAR and cannot be
    recovered. This moves a few multi-mapped reads (about 20 per sample on
    LJQQSK) and so changes about 0.1% of gene counts by a read or so. The test
    therefore requires:

      1. Portal, exact:     STAR uniquely and multi-mapped read counts. These
                            are computed before deduplication, so any difference
                            means the fastp/STAR recipe or reference is wrong.
      2. Portal, tolerance: dedup read counts and gene counts, within the
                            limits below. They are set well above the tie-order
                            noise but far below what a recipe error causes
                            (e.g. the wrong fastp settings left 18% of genes
                            differing).
      3. Repeat run, exact: with --repeat-out, a second run_pipeline.sh run on
                            the same input must be identical (STAR statistics,
                            dedup BAM records, counts). The pipeline writes
                            reads in FASTQ order, so results must not depend on
                            the run or the thread count.

    With --recipe-out, a validation/reproduce_portal.sh run is also compared,
    for information only. It used STAR's default output order, so it has its
    own tie-order noise and is not expected to match exactly.

Usage:
    python3 validation/compare_to_portal.py \
        --pipeline-out validation/work/det_a/out \
        --portal data/plasmidsaurus/LJQQSK \
        --run LJQQSK \
        [--repeat-out validation/work/det_b/out] \
        [--recipe-out validation/LJQQSK]

Inputs:
    --pipeline-out  run_pipeline.sh output directory (02_aligned/, 04_dedup/,
                    07_counts/gene_counts.txt)
    --portal        Portal download: results/<RUN>-expression-matrix.tsv,
                    results/<RUN>-mapping-stats-reads.csv, <RUN>_bam/*.bam
    --run           Plasmidsaurus run ID (file name prefix)
    --repeat-out    Optional second run_pipeline.sh output on the same input
    --recipe-out    Optional reproduce_portal.sh output (counts.txt)

Outputs:
    Report on stdout. Exit status 0 if every check passes, 1 otherwise.

Dependencies:
    Python >= 3.8 with pandas, numpy; samtools on PATH (use `pixi run`)

Author: Brent Chick
Date: 2026-09-26
Version: 2.0.0
"""

import argparse
import glob
import hashlib
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd


# ==============================================================================
# TOLERANCES VS THE PORTAL
# ==============================================================================
# Each limit is several times the worst value seen on LJQQSK (2026-09-26,
# 4 samples, both run_pipeline.sh and reproduce_portal.sh), shown in brackets.

# Dedup unique and multi-mapped read counts, parts per million [17 ppm]
DEDUP_MAX_PPM = 100
# Genes whose count differs at all, percent of all genes [0.12%]
GENES_DIFFERING_MAX_PCT = 0.5
# Largest difference for any single gene, reads [5]
GENE_MAX_ABS_DIFF = 10
# Total assigned counts, parts per million [1.0 ppm]
TOTAL_COUNTS_MAX_PPM = 10
# Pearson correlation of log1p counts [0.9999998]
PEARSON_MIN = 0.99999

# Counts are fractional (--fraction); treat differences below this as equal
EPS = 1e-6


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


def bam_records_md5(bam):
    """MD5 of a BAM's alignment records in file order, excluding the header.

    The header holds @PG command lines with run-specific paths, so it is left
    out; identical records in identical order means an identical result.
    """
    digest = hashlib.md5()
    with subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE) as proc:
        for chunk in iter(lambda: proc.stdout.read(1 << 20), b""):
            digest.update(chunk)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools view failed on {bam}")
    return digest.hexdigest()


def star_stats(out_dir, sample_id):
    """(uniquely mapped, multi-mapped) read counts from STAR's Log.final.out."""
    text = open(f"{out_dir}/02_aligned/{sample_id}_Log.final.out").read()
    get = lambda key: int(re.search(rf"{key} \|\s+(\d+)", text).group(1))
    return get("Uniquely mapped reads number"), get("Number of reads mapped to multiple loci")


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


def ppm(a, b):
    """Difference of a from b in parts per million of b."""
    return (a - b) / b * 1e6 if b else float("inf")


class Checks:
    """Collects pass/fail results and formats the status column."""

    def __init__(self):
        self.failures = []

    def __call__(self, ok, label):
        if not ok:
            self.failures.append(label)
        return "PASS" if ok else "FAIL"


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pipeline-out", required=True)
    ap.add_argument("--portal", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--repeat-out")
    ap.add_argument("--recipe-out")
    args = ap.parse_args()

    out, portal, run = args.pipeline_out, args.portal, args.run
    stats = pd.read_csv(f"{portal}/results/{run}-mapping-stats-reads.csv", index_col=0)
    ref = pd.read_csv(f"{portal}/results/{run}-expression-matrix.tsv", sep="\t", index_col=0)

    # The portal matrix lists every annotated gene; reindex so a gene missing
    # from ours shows up as a count difference rather than a crash
    def load_counts(path):
        return read_featurecounts(path).reindex(ref.index).fillna(0)

    ours = load_counts(f"{out}/07_counts/gene_counts.txt")
    repeat = load_counts(f"{args.repeat_out}/07_counts/gene_counts.txt") if args.repeat_out else None
    recipe = load_counts(f"{args.recipe_out}/counts.txt") if args.recipe_out else None

    check = Checks()
    print(f"=== {run}: run_pipeline.sh output vs portal ===")
    print(f"pipeline output: {out}")
    if args.repeat_out:
        print(f"repeat run:      {args.repeat_out}")
    print(f"tolerances: dedup {DEDUP_MAX_PPM} ppm; genes differing {GENES_DIFFERING_MAX_PCT}%; "
          f"max gene diff {GENE_MAX_ABS_DIFF} reads; total {TOTAL_COUNTS_MAX_PPM} ppm; "
          f"Pearson >= {PEARSON_MIN}")

    for sid in sorted(ours.columns):
        short = run_prefix(sid)
        # Portal mapping-stats rows are labelled without the run prefix
        s = stats.loc[sid[len(short) + 1:]]
        uniq, multi = star_stats(out, sid)
        dedup_bam = f"{out}/04_dedup/{sid}.dedup.bam"
        pu, pm = primary_nh_counts(glob.glob(f"{portal}/{run}_bam/{sid}_dedup-mapped-reads.bam")[0])
        ou, om = primary_nh_counts(dedup_bam)

        c, r = ours[sid], ref[f"{short}_count"]
        diff = (c - r).abs()
        n_diff = int((diff > EPS).sum())
        pct_diff = n_diff / len(diff) * 100
        pearson = np.corrcoef(np.log1p(c), np.log1p(r))[0, 1]
        total_ppm = ppm(c.sum(), r.sum())

        print(f"\n{sid}")
        print("  vs portal (exact)")
        print(f"    STAR unique    ours {uniq:>10}  portal {s['Uniquely Mapped']:>10}  "
              f"{uniq - s['Uniquely Mapped']:+d}  {check(uniq == s['Uniquely Mapped'], f'{sid}: STAR unique')}")
        print(f"    STAR multi     ours {multi:>10}  portal {s['Multi-mapped']:>10}  "
              f"{multi - s['Multi-mapped']:+d}  {check(multi == s['Multi-mapped'], f'{sid}: STAR multi')}")
        print("  vs portal (tolerance)")
        print(f"    dedup unique   ours {ou:>10}  portal {pu:>10}  {ou - pu:+d} ({ppm(ou, pu):+.1f} ppm)  "
              f"{check(abs(ppm(ou, pu)) <= DEDUP_MAX_PPM, f'{sid}: dedup unique')}")
        print(f"    dedup multi    ours {om:>10}  portal {pm:>10}  {om - pm:+d} ({ppm(om, pm):+.1f} ppm)  "
              f"{check(abs(ppm(om, pm)) <= DEDUP_MAX_PPM, f'{sid}: dedup multi')}")
        print(f"    counts total   ours {c.sum():>10.1f}  portal {r.sum():>10.1f}  ({total_ppm:+.2f} ppm)  "
              f"{check(abs(total_ppm) <= TOTAL_COUNTS_MAX_PPM, f'{sid}: total counts')}")
        print(f"    genes differing {n_diff} of {len(diff)} ({pct_diff:.3f}%)  "
              f"{check(pct_diff <= GENES_DIFFERING_MAX_PCT, f'{sid}: genes differing')}")
        print(f"    max gene diff  {diff.max():.3f} reads  "
              f"{check(diff.max() <= GENE_MAX_ABS_DIFF, f'{sid}: max gene diff')}")
        print(f"    Pearson r (log1p) {pearson:.7f}  {check(pearson >= PEARSON_MIN, f'{sid}: Pearson')}")

        if repeat is not None:
            ru, rm = star_stats(args.repeat_out, sid)
            rmax = (c - repeat[sid]).abs().max()
            same_bam = bam_records_md5(dedup_bam) == bam_records_md5(
                f"{args.repeat_out}/04_dedup/{sid}.dedup.bam")
            print("  vs repeat run (exact)")
            print(f"    STAR stats     {'identical' if (ru, rm) == (uniq, multi) else 'DIFFERENT'}  "
                  f"{check((ru, rm) == (uniq, multi), f'{sid}: repeat STAR stats')}")
            print(f"    dedup BAM      {'identical records' if same_bam else 'DIFFERENT records'}  "
                  f"{check(same_bam, f'{sid}: repeat dedup BAM')}")
            print(f"    counts         max |diff| {rmax:.3g}  "
                  f"{check(rmax <= EPS, f'{sid}: repeat counts')}")

        if recipe is not None:
            d = (recipe[short] - r).abs()
            print("  reproduce_portal.sh run vs portal (information only)")
            print(f"    genes differing {int((d > EPS).sum())}  max gene diff {d.max():.3f} reads  "
                  f"total {ppm(recipe[short].sum(), r.sum()):+.2f} ppm")

    print()
    if check.failures:
        print(f"RESULT: FAIL ({len(check.failures)} checks)")
        for f in check.failures:
            print(f"  - {f}")
        return 1
    print("RESULT: PASS (all checks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
