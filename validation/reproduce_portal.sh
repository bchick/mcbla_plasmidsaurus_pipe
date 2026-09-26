#!/usr/bin/env bash
#
# reproduce_portal.sh - Reproduce a Plasmidsaurus portal run from raw FASTQ
#
# Description:
#   Runs raw FASTQ -> fastp -> STAR -> sort -> UMICollapse -> featureCounts
#   with the settings recovered from Plasmidsaurus's own outputs (BAM @PG
#   headers, count matrix), and compares every stage against the portal files.
#   This is the concordance test for the pipeline: once it passes, the same
#   settings are carried into scripts/.
#
#   Each stage writes a .done marker, so re-running resumes where it stopped.
#
# Usage:
#   pixi run bash validation/reproduce_portal.sh
#
# Inputs (under data/):
#   plasmidsaurus/LJQQSK/fastq/*.fastq.gz          raw reads from the portal
#   plasmidsaurus/LJQQSK/LJQQSK_bam/*dedup*.bam    portal BAMs (header + stats)
#   plasmidsaurus/LJQQSK/results/                  portal count matrix + stats
#   reference/                                     Ensembl 114 FASTA/GTF, ERCC92
#
# Outputs:
#   validation/LJQQSK/                             all intermediate files
#   validation/LJQQSK/REPORT.txt                   stage-by-stage comparison
#
# Dependencies:
#   pixi environment (pixi.toml): STAR 2.7.11b, fastp 0.24.0, samtools 1.22.1,
#   UMICollapse 1.1.0, subread 2.1.1, python + pandas
#
# Author: Brent Chick
# Date: 2026-09-23

set -euo pipefail
[[ "${DEBUG:-}" == "true" ]] && set -x

# ==============================================================================
# CONFIGURATION
# ==============================================================================
readonly ROOT="$(cd "$(dirname "$0")/.." && pwd)"
readonly RUN="LJQQSK"
readonly PORTAL="${ROOT}/data/plasmidsaurus/${RUN}"
readonly REF="${ROOT}/data/reference"
readonly OUT="${ROOT}/validation/${RUN}"
readonly THREADS=32
readonly GENOME_FA="${REF}/Mus_musculus.GRCm39.114_ERCC.fa"
readonly GTF="${REF}/Mus_musculus.GRCm39.114_ERCC.gtf"
readonly STAR_INDEX="${REF}/Mus_musculus.GRCm39.114_ERCC.STAR"
readonly REPORT="${OUT}/REPORT.txt"

# Reads entering STAR for sample 1, from the portal's mapping-stats table.
# fastp settings are chosen by matching this number exactly.
readonly TARGET_STAR_INPUT=19715741

# Refuse to run outside the pixi environment. System copies of fastp/STAR
# (e.g. /usr/bin/fastp 0.20.1, STAR 2.7.11a) are on PATH on some hosts and
# would silently produce results that are not comparable to the portal.
[[ -n "${PIXI_PROJECT_ROOT:-}" ]] || {
    echo "ERROR: run inside the pixi environment: pixi run bash validation/reproduce_portal.sh" >&2
    exit 1
}
for tool in fastp STAR samtools umicollapse featureCounts; do
    command -v "${tool}" > /dev/null || { echo "ERROR: ${tool} not found in pixi env" >&2; exit 1; }
done

mkdir -p "${OUT}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
report() { echo "$*" | tee -a "${REPORT}" >&2; }
die() { log "ERROR: $*"; report "FAILED: $*"; exit 1; }
done_marker() { [[ -f "${OUT}/.$1.done" ]]; }
mark_done() { touch "${OUT}/.$1.done"; }

trap 'die "unexpected failure at line ${LINENO}"' ERR

report "=== ${RUN} portal reproduction, started $(date) ==="

# ==============================================================================
# STAGE 1: REFERENCE (Ensembl 114 primary assembly + ERCC92)
# ==============================================================================
# Plasmidsaurus's index is mus_musculus_GRCm39_114_ERCC. Check that the
# sequences we build from match their BAM header in name, length and order.
if ! done_marker reference; then
    log "Building genome FASTA"
    { zcat "${REF}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
        | sed -E 's/^>([^ ]+).*/>\1/'; cat "${REF}/ERCC92.fa"; } > "${GENOME_FA}"
    samtools faidx "${GENOME_FA}"

    portal_bam=$(ls "${PORTAL}"/LJQQSK_bam/*_dedup-mapped-reads.bam | head -1)
    samtools view -H "${portal_bam}" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $2"\t"$3}' \
        > "${OUT}/portal_sq.tsv"
    cut -f1,2 "${GENOME_FA}.fai" > "${OUT}/our_sq.tsv"
    if cmp -s "${OUT}/portal_sq.tsv" "${OUT}/our_sq.tsv"; then
        report "Reference: MATCH ($(wc -l < "${OUT}/our_sq.tsv") sequences, same names/lengths/order)"
    elif diff <(sort "${OUT}/portal_sq.tsv") <(sort "${OUT}/our_sq.tsv") > /dev/null; then
        report "Reference: same sequences, DIFFERENT ORDER (does not affect counts)"
    else
        report "Reference: MISMATCH vs portal BAM header:"
        diff "${OUT}/portal_sq.tsv" "${OUT}/our_sq.tsv" | head -20 | tee -a "${REPORT}" >&2 || true
    fi
    mark_done reference
fi

# ==============================================================================
# STAGE 2: CHOOSE FASTP SETTINGS
# ==============================================================================
# The portal only documents fastp loosely (poly-X trim, 3' quality trim, Q15,
# min length 50). Try the plausible readings on sample 1 and keep the one whose
# output read count equals the number of reads STAR saw.
readonly S1_FASTQ=$(ls "${PORTAL}"/fastq/${RUN}_1_*.fastq.gz)
declare -A FASTP_VARIANTS=(
    [A_polyx_cuttail_q15_len50]="--trim_poly_x --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 15 --qualified_quality_phred 15 --length_required 50"
    [B_polyx_q15_len50]="--trim_poly_x --qualified_quality_phred 15 --length_required 50"
    [C_polyx_cuttail_len50]="--trim_poly_x --cut_tail --cut_tail_mean_quality 15 --length_required 50"
    # E: fastp defaults for --cut_tail (window 4, mean Q20). Verified 2026-09-24 to
    # reproduce the portal read count exactly and per-read trimmed lengths 100%.
    [E_polyx_cuttail_default_q15_len50]="--trim_poly_x --cut_tail --qualified_quality_phred 15 --length_required 50"
    [D_polyx_cuttail_q15_len50_polyg]="--trim_poly_g --trim_poly_x --cut_tail --cut_tail_mean_quality 15 --qualified_quality_phred 15 --length_required 50"
)
if ! done_marker fastp_choice; then
    mkdir -p "${OUT}/fastp_variants"
    report ""
    report "fastp variants on sample 1 (target ${TARGET_STAR_INPUT} reads out):"
    best=""; best_diff=""
    for v in "${!FASTP_VARIANTS[@]}"; do
        # shellcheck disable=SC2086
        fastp -i "${S1_FASTQ}" -o /dev/null -w 16 ${FASTP_VARIANTS[$v]} \
            -j "${OUT}/fastp_variants/${v}.json" -h /dev/null 2> "${OUT}/fastp_variants/${v}.log"
        n=$(python3 -c "import json;print(json.load(open('${OUT}/fastp_variants/${v}.json'))['summary']['after_filtering']['total_reads'])")
        d=$(( n > TARGET_STAR_INPUT ? n - TARGET_STAR_INPUT : TARGET_STAR_INPUT - n ))
        report "  ${v}: ${n} (off by ${d})"
        if [[ -z "${best}" ]] || (( d < best_diff )); then best="${v}"; best_diff="${d}"; fi
    done
    echo "${best}" > "${OUT}/fastp_choice.txt"
    report "  -> using ${best}$([[ ${best_diff} -eq 0 ]] && echo ' (EXACT)' || echo " (closest; off by ${best_diff})")"
    mark_done fastp_choice
fi
readonly FASTP_ARGS="${FASTP_VARIANTS[$(cat "${OUT}/fastp_choice.txt")]}"

# ==============================================================================
# STAGE 3: STAR INDEX
# ==============================================================================
# Defaults except for the annotation; --sjdbOverhang 100 is STAR's default.
if ! done_marker star_index; then
    log "Building STAR index (about an hour)"
    mkdir -p "${STAR_INDEX}"
    STAR --runMode genomeGenerate --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GTF}" \
        --outFileNamePrefix "${STAR_INDEX}/" > "${OUT}/star_index.log" 2>&1
    mark_done star_index
fi

# ==============================================================================
# STAGE 4: PER-SAMPLE PROCESSING
# ==============================================================================
for fq in "${PORTAL}"/fastq/*.fastq.gz; do
    name=$(basename "${fq}" .fastq.gz)          # e.g. LJQQSK_1_BR4_1_US_r1
    sid=$(cut -d_ -f1,2 <<< "${name}")          # e.g. LJQQSK_1
    d="${OUT}/${sid}"; mkdir -p "${d}"

    if ! done_marker "${sid}.fastp"; then
        log "${sid}: fastp"
        # shellcheck disable=SC2086
        fastp -i "${fq}" -o "${d}/trimmed.fastq.gz" -w 16 ${FASTP_ARGS} \
            -j "${d}/fastp.json" -h "${d}/fastp.html" 2> "${d}/fastp.log"
        mark_done "${sid}.fastp"
    fi

    if ! done_marker "${sid}.star"; then
        log "${sid}: STAR"
        # Exactly the portal's STAR command (from the BAM @PG header), minus
        # --quantMode TranscriptomeSAM, which only adds a transcriptome BAM.
        STAR --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${d}/trimmed.fastq.gz" \
            --readFilesCommand pigz -dc -p 8 \
            --outFileNamePrefix "${d}/" \
            --outSAMtype BAM Unsorted \
            --outSAMattributes NH HI AS nM NM MD \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outReadsUnmapped Fastx > "${d}/star.stdout" 2>&1
        samtools sort -@ 16 -o "${d}/sorted.bam" "${d}/Aligned.out.bam"
        samtools index "${d}/sorted.bam"
        rm -f "${d}/Aligned.out.bam"
        mark_done "${sid}.star"
    fi

    if ! done_marker "${sid}.dedup"; then
        log "${sid}: UMICollapse"
        # The UMI follows the last "_" in the read name (e.g. ..._CACTTGCGCGAGCA)
        _JAVA_OPTIONS="-Xmx64g -Xss1g" umicollapse bam \
            -i "${d}/sorted.bam" -o "${d}/dedup.bam" --umi-sep _ \
            > "${d}/umicollapse.log" 2>&1
        samtools index "${d}/dedup.bam"
        mark_done "${sid}.dedup"
    fi
done

# ==============================================================================
# STAGE 5: COUNTING
# ==============================================================================
# Settings reproduced exactly from the portal BAMs: whole gene bodies, forward
# strand, multi-mappers split fractionally.
if ! done_marker counts; then
    log "featureCounts"
    featureCounts -T "${THREADS}" -a "${GTF}" -o "${OUT}/counts.txt" \
        -t gene -g gene_id -s 1 -M --fraction \
        "${OUT}"/${RUN}_*/dedup.bam > "${OUT}/featurecounts.log" 2>&1
    mark_done counts
fi

# ==============================================================================
# STAGE 6: COMPARISON REPORT
# ==============================================================================
python3 - "${OUT}" "${PORTAL}" "${RUN}" <<'PYEOF' 2>&1 | tee -a "${REPORT}"
import sys, re, glob, subprocess, pandas as pd, numpy as np
out, portal, run = sys.argv[1:]
stats = pd.read_csv(f"{portal}/results/{run}-mapping-stats-reads.csv", index_col=0)
ref = pd.read_csv(f"{portal}/results/{run}-expression-matrix.tsv", sep="\t", index_col=0)
ours = pd.read_csv(f"{out}/counts.txt", sep="\t", comment="#", index_col=0).iloc[:, 5:]
ours.columns = [c.split("/")[-2] for c in ours.columns]
ours = ours.reindex(ref.index).fillna(0)

def primary_nh_counts(bam):
    """Primary mapped reads split by NH==1 (unique) vs NH>1 (multi)."""
    p = subprocess.run(f"samtools view -F 260 {bam} | grep -o 'NH:i:[0-9]*' | sort | uniq -c",
                       shell=True, capture_output=True, text=True).stdout
    c = {int(nh.split(':')[-1]): int(n) for n, nh in (l.split() for l in p.strip().splitlines())}
    return c.get(1, 0), sum(v for k, v in c.items() if k > 1)

def pct(a, b): return f"{(a - b) / b * 100:+.2f}%"

print("\nPer-sample comparison (ours vs portal)")
for sid in ours.columns:
    log = open(f"{out}/{sid}/Log.final.out").read()
    get = lambda k: int(re.search(rf"{k} \|\s+(\d+)", log).group(1))
    inp, uniq, multi = get("Number of input reads"), get("Uniquely mapped reads number"), get("Number of reads mapped to multiple loci")
    label = next(i for i in stats.index if sid.split('_')[1] == i.split('_')[1])
    s = stats.loc[label]
    pbam = glob.glob(f"{portal}/{run}_bam/{sid}_*_dedup-mapped-reads.bam")[0]
    pu, pm = primary_nh_counts(pbam)
    ou, om = primary_nh_counts(f"{out}/{sid}/dedup.bam")
    c, r = ours[sid], ref[f"{sid}_count"]
    diff = (c - r).abs()
    print(f"\n{sid} ({label})")
    print(f"  STAR unique      ours {uniq:>10}  portal {s['Uniquely Mapped']:>10}  {pct(uniq, s['Uniquely Mapped'])}")
    print(f"  STAR multi       ours {multi:>10}  portal {s['Multi-mapped']:>10}  {pct(multi, s['Multi-mapped'])}")
    print(f"  dedup unique     ours {ou:>10}  portal {pu:>10}  {pct(ou, pu)}")
    print(f"  dedup multi      ours {om:>10}  portal {pm:>10}  {pct(om, pm)}")
    print(f"  counts total     ours {c.sum():>10.0f}  portal {r.sum():>10.0f}  {pct(c.sum(), r.sum())}")
    print(f"  genes exact {np.mean(diff < 1e-6)*100:.2f}%   within 1 read {np.mean(diff <= 1)*100:.2f}%   "
          f"Pearson r (log1p) {np.corrcoef(np.log1p(c), np.log1p(r))[0,1]:.5f}")
PYEOF

report ""
report "=== finished $(date) ==="
