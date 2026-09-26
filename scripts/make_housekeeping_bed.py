#!/usr/bin/env python3
"""
make_housekeeping_bed.py - Build a housekeeping-gene BED12 for RSeQC gene body coverage

Description:
    RSeQC's geneBody_coverage.py scales with the number of transcripts in its
    BED file. On the full Ensembl annotation (~278k transcripts) it takes hours
    per sample, so RSeQC's authors recommend running it on housekeeping genes
    only: they are expressed in every tissue, so their 5'->3' coverage profile
    reflects RNA integrity rather than biology.

    RSeQC ships housekeeping lists as BED files in an older assembly with RefSeq
    transcript IDs (e.g. mm10.HouseKeepingGenes.bed). Their coordinates do not
    match GRCm39, so this script keeps only the RefSeq IDs, maps them to Ensembl
    genes with Ensembl's own RefSeq cross-reference table for the same release,
    and writes each gene's Ensembl canonical transcript from the pipeline's
    BED12. The output therefore uses exactly the reference the reads were
    aligned to.

Usage:
    python3 scripts/make_housekeeping_bed.py \
        --housekeeping mm10.HouseKeepingGenes.bed.gz \
        --refseq-xref Mus_musculus.GRCm39.114.refseq.tsv.gz \
        --gtf data/reference/Mus_musculus.GRCm39.114_ERCC.gtf \
        --bed12 data/reference/Mus_musculus.GRCm39.114_ERCC.bed12 \
        --output data/reference/Mus_musculus.GRCm39.114.housekeeping.bed12

Inputs:
    --housekeeping  RSeQC housekeeping BED (column 4 = RefSeq transcript ID);
                    https://sourceforge.net/projects/rseqc/files/BED/
    --refseq-xref   Ensembl <species>.<assembly>.<release>.refseq.tsv(.gz), from
                    https://ftp.ensembl.org/pub/release-<N>/tsv/<species>/
    --gtf           Ensembl GTF used to build the STAR index (plain or .gz)
    --bed12         BED12 made from that GTF (column 4 = Ensembl transcript ID)

Outputs:
    --output        BED12 of one canonical transcript per housekeeping gene.
                    A summary of how many IDs mapped is written to stderr.

Dependencies:
    Python >= 3.8 (standard library only)

Author: Brent Chick
Date: 2026-09-26
Version: 1.0.0
"""

import argparse
import gzip
import re
import sys


# ==============================================================================
# HELPERS
# ==============================================================================

def open_text(path):
    """Open a plain or gzipped text file for reading."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def log(message):
    print(message, file=sys.stderr)


def read_housekeeping_ids(path):
    """RefSeq transcript IDs (version stripped) from column 4 of a BED file."""
    with open_text(path) as fh:
        return {line.split("\t")[3].split(".")[0]
                for line in fh if line.strip() and not line.startswith(("#", "track"))}


def refseq_to_gene(path, wanted):
    """Map RefSeq transcript IDs to Ensembl gene IDs.

    Only curated mRNA cross-references (NM_) are used; predicted models (XM_)
    are not in the RSeQC lists. One RefSeq ID can map to several Ensembl genes
    (e.g. duplicated loci); all are kept, since each is a housekeeping locus.
    """
    mapping = {}
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gene_col, xref_col = header.index("gene_stable_id"), header.index("xref")
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            xref = fields[xref_col].split(".")[0]
            if xref in wanted:
                mapping.setdefault(xref, set()).add(fields[gene_col])
    return mapping


def canonical_transcripts(gtf_path, genes):
    """Ensembl canonical transcript ID for each gene in `genes`.

    Ensembl tags exactly one transcript per gene as "Ensembl_canonical"; it is
    the most representative model (MANE Select where one exists), so using it
    avoids weighting genes by their number of isoforms.
    """
    gene_re = re.compile(r'gene_id "([^"]+)"')
    tx_re = re.compile(r'transcript_id "([^"]+)"')
    canonical = {}
    with open_text(gtf_path) as fh:
        for line in fh:
            # Cheap substring tests first: the GTF is over 1 GB
            if "Ensembl_canonical" not in line or "\ttranscript\t" not in line:
                continue
            gene = gene_re.search(line).group(1)
            if gene in genes:
                canonical[gene] = tx_re.search(line).group(1)
    return canonical


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--housekeeping", required=True)
    ap.add_argument("--refseq-xref", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--bed12", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    hk_ids = read_housekeeping_ids(args.housekeeping)
    log(f"Housekeeping RefSeq IDs: {len(hk_ids)}")

    xref = refseq_to_gene(args.refseq_xref, hk_ids)
    genes = set().union(*xref.values()) if xref else set()
    log(f"  mapped to Ensembl genes: {len(xref)} IDs -> {len(genes)} genes "
        f"({len(hk_ids) - len(xref)} IDs not in this Ensembl release)")

    canonical = canonical_transcripts(args.gtf, genes)
    wanted_tx = set(canonical.values())
    log(f"  genes with a canonical transcript: {len(canonical)}")

    written = 0
    with open_text(args.bed12) as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.split("\t", 4)[3] in wanted_tx:
                fout.write(line)
                written += 1
    log(f"Wrote {written} transcripts to {args.output}")

    # Far fewer than the ~3,800 genes in RSeQC's lists means the inputs do not
    # belong together (wrong species, or a BED12 not made from this GTF)
    if written < 1000:
        log("ERROR: fewer than 1000 transcripts written; check that the inputs "
            "are for the same species and Ensembl release")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
