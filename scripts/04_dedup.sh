#!/usr/bin/env bash
#
# 04_dedup.sh - Remove PCR and optical duplicates using UMICollapse
#
# Description:
#   This script performs UMI-based deduplication of aligned reads using UMICollapse.
#   UMI (Unique Molecular Identifier) deduplication is essential for accurate
#   quantification in RNA-seq, as it removes PCR duplicates that can bias expression
#   estimates. Unlike position-only deduplication, UMI-based methods use molecular
#   barcodes to distinguish true biological duplicates from PCR artifacts.
#
# Usage:
#   ./04_dedup.sh -i <input_bam> -o <output_dir> [-s sample_id]
#
# Arguments:
#   -i, --input         Path to coordinate-sorted BAM file (required)
#   -o, --output-dir    Output directory for deduplicated BAM (required)
#   -s, --sample-id     Sample identifier (default: derived from input filename)
#   -u, --umi-sep       UMI separator in read name (default: :)
#   -a, --algorithm     Deduplication algorithm (default: directional)
#   -h, --help          Display this help message
#
# Dependencies:
#   - UMICollapse >= 1.1.0 (UMI-based deduplication)
#   - samtools >= 1.22.1 (BAM indexing)
#
# Output:
#   - {sample}.dedup.bam       Deduplicated BAM file
#   - {sample}.dedup.bam.bai   BAM index file
#   - {sample}.dedup_stats.txt Deduplication statistics
#
# Example:
#   ./04_dedup.sh -i sample.sorted.bam -o results/04_dedup/
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
readonly DEFAULT_UMI_SEP=":"
readonly DEFAULT_ALGORITHM="directional"
readonly DEFAULT_EDIT_DISTANCE=1

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
Usage: ${SCRIPT_NAME} -i <input_bam> -o <output_dir> [OPTIONS]

Remove PCR and optical duplicates using UMI-based deduplication.

Required Arguments:
  -i, --input         Path to coordinate-sorted BAM file
  -o, --output-dir    Output directory for deduplicated BAM

Optional Arguments:
  -s, --sample-id     Sample identifier (default: derived from input filename)
  -u, --umi-sep       UMI separator in read name (default: ${DEFAULT_UMI_SEP})
  -a, --algorithm     Deduplication algorithm: directional, adjacency, cluster
                      (default: ${DEFAULT_ALGORITHM})
  -h, --help          Display this help message

Algorithms:
  directional  - Considers UMIs within edit distance and read direction
  adjacency    - Groups UMIs within edit distance threshold
  cluster      - Hierarchical clustering of UMIs

Examples:
  ${SCRIPT_NAME} -i sample.sorted.bam -o results/04_dedup/
  ${SCRIPT_NAME} -i sample.sorted.bam -o results/04_dedup/ -a adjacency

EOF
    exit 0
}

check_command() {
    local cmd="$1"
    if ! command -v "${cmd}" &> /dev/null; then
        die "Required command '${cmd}' not found in PATH"
    fi
}

extract_sample_id() {
    local filename="$1"
    filename="$(basename "${filename}")"
    filename="${filename%.bam}"
    filename="${filename%.sorted}"
    echo "${filename}"
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

input_bam=""
output_dir=""
sample_id=""
umi_sep="${DEFAULT_UMI_SEP}"
algorithm="${DEFAULT_ALGORITHM}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input_bam="$2"
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
        -u|--umi-sep)
            umi_sep="$2"
            shift 2
            ;;
        -a|--algorithm)
            algorithm="$2"
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
[[ -z "${output_dir}" ]] && die "Output directory is required (-o)"
[[ -f "${input_bam}" ]] || die "Input BAM file not found: ${input_bam}"
[[ -f "${input_bam}.bai" ]] || die "BAM index not found: ${input_bam}.bai (run samtools index first)"

if [[ -z "${sample_id}" ]]; then
    sample_id="$(extract_sample_id "${input_bam}")"
    log "INFO" "Derived sample ID: ${sample_id}"
fi

# Validate algorithm choice
case "${algorithm}" in
    directional|adjacency|cluster)
        ;;
    *)
        die "Invalid algorithm: ${algorithm}. Choose from: directional, adjacency, cluster"
        ;;
esac

check_command "java"
check_command "samtools"
mkdir -p "${output_dir}"

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Starting UMI deduplication for: ${sample_id}"
log "INFO" "=========================================="

log "INFO" "Algorithm: ${algorithm}"
log "INFO" "UMI separator: '${umi_sep}'"

output_bam="${output_dir}/${sample_id}.dedup.bam"
stats_file="${output_dir}/${sample_id}.dedup_stats.txt"

# Count reads before deduplication for statistics
reads_before=$(samtools view -c "${input_bam}")
log "INFO" "Reads before deduplication: ${reads_before}"

# ==============================================================================
# UMI DEDUPLICATION
# ==============================================================================
log "INFO" "Running UMICollapse..."

# UMICollapse parameters:
# bam: Input/output mode for BAM files
# -i: Input BAM file
# -o: Output BAM file
# --umi-sep: Character separating UMI from read name
# --algo: Algorithm for UMI grouping
# --edit-distance: Maximum edit distance for UMI comparison
#
# The directional algorithm is recommended for most RNA-seq applications as it
# accounts for the expected error profile in UMI sequences and read direction

# Check if UMICollapse is available as a JAR or command
if command -v umicollapse &> /dev/null; then
    # UMICollapse available as command
    umicollapse bam \
        -i "${input_bam}" \
        -o "${output_bam}" \
        --umi-sep "${umi_sep}" \
        --algo "${algorithm}" \
        --edit-distance "${DEFAULT_EDIT_DISTANCE}" \
        2>&1 | tee "${stats_file}"
elif [[ -n "${UMICOLLAPSE_JAR:-}" ]] && [[ -f "${UMICOLLAPSE_JAR}" ]]; then
    # UMICollapse available as JAR file
    java -jar "${UMICOLLAPSE_JAR}" bam \
        -i "${input_bam}" \
        -o "${output_bam}" \
        --umi-sep "${umi_sep}" \
        --algo "${algorithm}" \
        --edit-distance "${DEFAULT_EDIT_DISTANCE}" \
        2>&1 | tee "${stats_file}"
else
    # Fallback: try to find UMICollapse wrapper script or JAR in common locations
    umicollapse_script=""
    for script_path in \
        "/data/bchick/tools/umicollapse/umicollapse.sh" \
        "${HOME}/tools/umicollapse/umicollapse.sh" \
        "/opt/umicollapse/umicollapse.sh"; do
        if [[ -f "${script_path}" ]]; then
            umicollapse_script="${script_path}"
            break
        fi
    done

    if [[ -n "${umicollapse_script}" ]]; then
        log "INFO" "Found UMICollapse: ${umicollapse_script}"
        "${umicollapse_script}" bam \
            -i "${input_bam}" \
            -o "${output_bam}" \
            --umi-sep "${umi_sep}" \
            --algo "${algorithm}" \
            --edit-distance "${DEFAULT_EDIT_DISTANCE}" \
            2>&1 | tee "${stats_file}"
    else
        die "UMICollapse not found. Install to /data/bchick/tools/umicollapse/ or set UMICOLLAPSE_JAR environment variable."
    fi
fi

log "INFO" "UMICollapse completed"

# ==============================================================================
# BAM INDEXING
# ==============================================================================
log "INFO" "Creating BAM index..."

samtools index "${output_bam}"

# ==============================================================================
# OUTPUT VERIFICATION AND STATISTICS
# ==============================================================================
log "INFO" "Verifying outputs..."

[[ -s "${output_bam}" ]] || die "Output BAM file is empty or missing: ${output_bam}"
[[ -s "${output_bam}.bai" ]] || die "BAM index file is empty or missing: ${output_bam}.bai"

# Count reads after deduplication
reads_after=$(samtools view -c "${output_bam}")

# Calculate deduplication rate
if [[ "${reads_before}" -gt 0 ]]; then
    duplicates_removed=$((reads_before - reads_after))
    dedup_rate=$(echo "scale=2; ${duplicates_removed} * 100 / ${reads_before}" | bc)

    log "INFO" "Deduplication statistics:"
    log "INFO" "  Reads before:        ${reads_before}"
    log "INFO" "  Reads after:         ${reads_after}"
    log "INFO" "  Duplicates removed:  ${duplicates_removed}"
    log "INFO" "  Deduplication rate:  ${dedup_rate}%"

    # Append summary to stats file
    {
        echo ""
        echo "=== Summary ==="
        echo "Reads before deduplication: ${reads_before}"
        echo "Reads after deduplication: ${reads_after}"
        echo "Duplicates removed: ${duplicates_removed}"
        echo "Deduplication rate: ${dedup_rate}%"
    } >> "${stats_file}"
fi

# ==============================================================================
# COMPLETION
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "UMI deduplication completed for: ${sample_id}"
log "INFO" "Output BAM: ${output_bam}"
log "INFO" "=========================================="

exit 0
