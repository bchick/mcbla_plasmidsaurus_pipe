#!/usr/bin/env bash
#
# 05_mapping_qc.sh - Generate mapping quality metrics using RSeQC and Qualimap
#
# Description:
#   This script generates comprehensive alignment quality metrics using RSeQC
#   and Qualimap. These tools assess alignment quality, strand specificity,
#   read distribution across genomic features, and other key QC metrics essential
#   for evaluating RNA-seq data quality.
#
# Usage:
#   ./05_mapping_qc.sh -i <input_bam> -b <bed_file> -o <output_dir> [-g gtf_file]
#
# Arguments:
#   -i, --input         Path to deduplicated BAM file (required)
#   -b, --bed           Path to gene model BED file for RSeQC (required)
#   -g, --gtf           Path to GTF annotation file for Qualimap (optional)
#   -o, --output-dir    Output directory for QC reports (required)
#   -s, --sample-id     Sample identifier (default: derived from input filename)
#   -t, --threads       Number of threads (default: 8)
#   -h, --help          Display this help message
#
# Dependencies:
#   - RSeQC >= 5.0.4 (RNA-seq quality control)
#   - Qualimap >= 2.3 (general alignment QC)
#
# Output:
#   - {sample}/rseqc/           RSeQC output directory
#   - {sample}/qualimap/        Qualimap output directory
#
# Example:
#   ./05_mapping_qc.sh -i sample.dedup.bam -b genes.bed -g genes.gtf \
#       -o results/05_qc/ -t 8
#
# Author: Plasmidsaurus RNA-seq Pipeline
# Date: 2024
# Version: 1.0.0

# ==============================================================================
# SHELL OPTIONS
# ==============================================================================
set -euo pipefail
[[ "${DEBUG:-}" == "true" ]] && set -x

# ==============================================================================
# CONSTANTS
# ==============================================================================
readonly SCRIPT_NAME="$(basename "$0")"
readonly DEFAULT_THREADS=8

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

log() {
    local level="$1"
    shift
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level}] $*" >&2
}

die() {
    log "ERROR" "$*"
    exit 1
}

usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} -i <input_bam> -b <bed_file> -o <output_dir> [OPTIONS]

Generate mapping quality metrics using RSeQC and Qualimap.

Required Arguments:
  -i, --input         Path to deduplicated BAM file
  -b, --bed           Path to gene model BED file for RSeQC
  -o, --output-dir    Output directory for QC reports

Optional Arguments:
  -g, --gtf           Path to GTF annotation file for Qualimap
  -s, --sample-id     Sample identifier (default: derived from input filename)
  -t, --threads       Number of threads (default: ${DEFAULT_THREADS})
  -h, --help          Display this help message

Examples:
  ${SCRIPT_NAME} -i sample.dedup.bam -b genes.bed -o results/05_qc/
  ${SCRIPT_NAME} -i sample.dedup.bam -b genes.bed -g genes.gtf -o results/05_qc/ -t 16

EOF
    exit 0
}

check_command() {
    local cmd="$1"
    if ! command -v "${cmd}" &> /dev/null; then
        log "WARN" "Command '${cmd}' not found - skipping related analyses"
        return 1
    fi
    return 0
}

extract_sample_id() {
    local filename="$1"
    filename="$(basename "${filename}")"
    filename="${filename%.bam}"
    filename="${filename%.dedup}"
    filename="${filename%.sorted}"
    echo "${filename}"
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

input_bam=""
bed_file=""
gtf_file=""
output_dir=""
sample_id=""
threads="${DEFAULT_THREADS}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input_bam="$2"
            shift 2
            ;;
        -b|--bed)
            bed_file="$2"
            shift 2
            ;;
        -g|--gtf)
            gtf_file="$2"
            shift 2
            ;;
        -o|--output-dir)
            output_dir="$2"
            shift 2
            ;;
        -s|--sample-id)
            sample_id="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            die "Unknown option: $1. Use -h for help."
            ;;
    esac
done

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================
log "INFO" "Validating inputs..."

[[ -z "${input_bam}" ]] && die "Input BAM file is required (-i)"
[[ -z "${bed_file}" ]] && die "BED file is required (-b)"
[[ -z "${output_dir}" ]] && die "Output directory is required (-o)"
[[ -f "${input_bam}" ]] || die "Input BAM file not found: ${input_bam}"
[[ -f "${bed_file}" ]] || die "BED file not found: ${bed_file}"

if [[ -n "${gtf_file}" ]] && [[ ! -f "${gtf_file}" ]]; then
    die "GTF file not found: ${gtf_file}"
fi

if [[ -z "${sample_id}" ]]; then
    sample_id="$(extract_sample_id "${input_bam}")"
    log "INFO" "Derived sample ID: ${sample_id}"
fi

mkdir -p "${output_dir}"

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Starting mapping QC for: ${sample_id}"
log "INFO" "=========================================="

# Create sample-specific output directories
rseqc_dir="${output_dir}/${sample_id}/rseqc"
qualimap_dir="${output_dir}/${sample_id}/qualimap"
mkdir -p "${rseqc_dir}" "${qualimap_dir}"

# ==============================================================================
# RSEQC QUALITY CONTROL
# ==============================================================================
# RSeQC provides RNA-seq specific QC metrics including:
# - Strand specificity (infer_experiment.py)
# - Read distribution across genomic features (read_distribution.py)
# - Gene body coverage (geneBody_coverage.py)
# - Inner distance (for paired-end, inner_distance.py)

if check_command "infer_experiment.py"; then
    log "INFO" "Running RSeQC analyses..."

    # ------------------------------------------------------------------
    # Infer Experiment: Determine strand specificity
    # ------------------------------------------------------------------
    # This is crucial for proper featureCounts configuration
    # Outputs fraction of reads mapping to sense/antisense strands
    log "INFO" "Running infer_experiment.py (strand specificity)..."
    infer_experiment.py \
        -r "${bed_file}" \
        -i "${input_bam}" \
        > "${rseqc_dir}/infer_experiment.txt" \
        2>> "${rseqc_dir}/rseqc.log"

    # ------------------------------------------------------------------
    # Read Distribution: Distribution across genomic features
    # ------------------------------------------------------------------
    # Shows percentage of reads mapping to CDS, UTRs, introns, intergenic
    log "INFO" "Running read_distribution.py (feature distribution)..."
    read_distribution.py \
        -r "${bed_file}" \
        -i "${input_bam}" \
        > "${rseqc_dir}/read_distribution.txt" \
        2>> "${rseqc_dir}/rseqc.log"

    # ------------------------------------------------------------------
    # BAM Statistics: Basic alignment statistics
    # ------------------------------------------------------------------
    log "INFO" "Running bam_stat.py (alignment statistics)..."
    bam_stat.py \
        -i "${input_bam}" \
        > "${rseqc_dir}/bam_stat.txt" \
        2>> "${rseqc_dir}/rseqc.log"

    # ------------------------------------------------------------------
    # Junction Saturation: Assess sequencing depth for splice junctions
    # ------------------------------------------------------------------
    log "INFO" "Running junction_saturation.py (junction saturation)..."
    junction_saturation.py \
        -r "${bed_file}" \
        -i "${input_bam}" \
        -o "${rseqc_dir}/junction_saturation" \
        2>> "${rseqc_dir}/rseqc.log" || \
        log "WARN" "junction_saturation.py failed - continuing"

    # ------------------------------------------------------------------
    # Gene Body Coverage: 5' to 3' coverage bias
    # ------------------------------------------------------------------
    # Important for detecting RNA degradation
    log "INFO" "Running geneBody_coverage.py (coverage bias)..."
    geneBody_coverage.py \
        -r "${bed_file}" \
        -i "${input_bam}" \
        -o "${rseqc_dir}/geneBody_coverage" \
        2>> "${rseqc_dir}/rseqc.log" || \
        log "WARN" "geneBody_coverage.py failed - continuing"

    log "INFO" "RSeQC analyses completed"
else
    log "WARN" "RSeQC not found - skipping RSeQC analyses"
fi

# ==============================================================================
# QUALIMAP QUALITY CONTROL
# ==============================================================================
# Qualimap provides general alignment QC metrics including:
# - Mapping quality distribution
# - Coverage statistics
# - GC content
# - Insert size distribution (for paired-end)

if check_command "qualimap"; then
    log "INFO" "Running Qualimap RNA-seq analysis..."

    # Build qualimap command
    qualimap_cmd=(
        qualimap rnaseq
        -bam "${input_bam}"
        -outdir "${qualimap_dir}"
        -outformat HTML
        --java-mem-size=8G
    )

    # Add GTF if provided (enables feature counting by Qualimap)
    if [[ -n "${gtf_file}" ]]; then
        qualimap_cmd+=(-gtf "${gtf_file}")
    fi

    # Execute Qualimap
    "${qualimap_cmd[@]}" \
        2>&1 | tee "${qualimap_dir}/qualimap.log"

    log "INFO" "Qualimap analysis completed"
else
    log "WARN" "Qualimap not found - skipping Qualimap analysis"
fi

# ==============================================================================
# QC SUMMARY
# ==============================================================================
log "INFO" "Generating QC summary..."

summary_file="${output_dir}/${sample_id}/qc_summary.txt"
{
    echo "=========================================="
    echo "QC Summary for: ${sample_id}"
    echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================="
    echo ""

    # Include strand specificity results if available
    if [[ -f "${rseqc_dir}/infer_experiment.txt" ]]; then
        echo "=== Strand Specificity ==="
        cat "${rseqc_dir}/infer_experiment.txt"
        echo ""
    fi

    # Include read distribution if available
    if [[ -f "${rseqc_dir}/read_distribution.txt" ]]; then
        echo "=== Read Distribution ==="
        cat "${rseqc_dir}/read_distribution.txt"
        echo ""
    fi

    # Include BAM statistics if available
    if [[ -f "${rseqc_dir}/bam_stat.txt" ]]; then
        echo "=== BAM Statistics ==="
        cat "${rseqc_dir}/bam_stat.txt"
        echo ""
    fi

} > "${summary_file}"

log "INFO" "QC summary written to: ${summary_file}"

# ==============================================================================
# COMPLETION
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Mapping QC completed for: ${sample_id}"
log "INFO" "RSeQC output: ${rseqc_dir}"
log "INFO" "Qualimap output: ${qualimap_dir}"
log "INFO" "=========================================="

exit 0
