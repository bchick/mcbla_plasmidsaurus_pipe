#!/usr/bin/env bash
#
# 03_sort_bam.sh - Coordinate sort BAM files using samtools
#
# Description:
#   This script performs coordinate sorting of BAM files using samtools. While
#   STAR can output coordinate-sorted BAMs directly, this step ensures proper
#   sorting and adds BAM indexing for downstream tools. The script also validates
#   the BAM file integrity.
#
#   Note: If STAR was run with --outSAMtype BAM SortedByCoordinate, this step
#   primarily adds the BAM index (.bai) and validates the file.
#
# Usage:
#   ./03_sort_bam.sh -i <input_bam> -o <output_dir> [-t threads] [-m memory]
#
# Arguments:
#   -i, --input         Path to input BAM file (required)
#   -o, --output-dir    Output directory for sorted BAM (required)
#   -s, --sample-id     Sample identifier (default: derived from input filename)
#   -t, --threads       Number of threads (default: 8)
#   -m, --memory        Memory per thread for sorting (default: 2G)
#   -h, --help          Display this help message
#
# Dependencies:
#   - samtools >= 1.22.1 (BAM manipulation)
#
# Output:
#   - {sample}.sorted.bam      Coordinate-sorted BAM file
#   - {sample}.sorted.bam.bai  BAM index file
#
# Example:
#   ./03_sort_bam.sh -i sample_Aligned.sortedByCoord.out.bam \
#       -o results/03_sorted/ -t 8 -m 4G
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
readonly DEFAULT_MEMORY="2G"

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

Coordinate sort BAM files and create index using samtools.

Required Arguments:
  -i, --input         Path to input BAM file
  -o, --output-dir    Output directory for sorted BAM

Optional Arguments:
  -s, --sample-id     Sample identifier (default: derived from input filename)
  -t, --threads       Number of threads (default: ${DEFAULT_THREADS})
  -m, --memory        Memory per thread for sorting (default: ${DEFAULT_MEMORY})
  -h, --help          Display this help message

Examples:
  ${SCRIPT_NAME} -i sample_Aligned.out.bam -o results/03_sorted/ -t 8 -m 4G

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
    filename="${filename%_Aligned.sortedByCoord.out}"
    filename="${filename%_Aligned.out}"
    filename="${filename%.sorted}"
    echo "${filename}"
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

input_bam=""
output_dir=""
sample_id=""
threads="${DEFAULT_THREADS}"
memory="${DEFAULT_MEMORY}"

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
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -m|--memory)
            memory="$2"
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

if [[ -z "${sample_id}" ]]; then
    sample_id="$(extract_sample_id "${input_bam}")"
    log "INFO" "Derived sample ID: ${sample_id}"
fi

check_command "samtools"
mkdir -p "${output_dir}"

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Starting BAM processing for: ${sample_id}"
log "INFO" "=========================================="

# Log tool version
samtools_version="$(samtools --version | head -1)"
log "INFO" "Using ${samtools_version}"

output_bam="${output_dir}/${sample_id}.sorted.bam"

# Check if input BAM is already coordinate sorted
# This avoids unnecessary re-sorting of STAR output
input_sort_order=$(samtools view -H "${input_bam}" | grep "^@HD" | grep -o "SO:[a-z]*" || echo "SO:unknown")
log "INFO" "Input BAM sort order: ${input_sort_order}"

if [[ "${input_sort_order}" == "SO:coordinate" ]]; then
    # BAM is already sorted, just copy and index
    log "INFO" "Input BAM is already coordinate sorted"
    log "INFO" "Copying BAM and creating index..."

    # Use samtools view to copy (validates the file in the process)
    samtools view \
        -@ "${threads}" \
        -b \
        -o "${output_bam}" \
        "${input_bam}"
else
    # Need to sort the BAM file
    log "INFO" "Sorting BAM by coordinate..."

    # samtools sort parameters:
    # -@ threads: Number of additional threads for compression
    # -m memory: Maximum memory per thread
    # -o output: Output file path
    samtools sort \
        -@ "${threads}" \
        -m "${memory}" \
        -o "${output_bam}" \
        "${input_bam}"
fi

log "INFO" "BAM sorting/copying completed"

# ==============================================================================
# BAM INDEXING
# ==============================================================================
log "INFO" "Creating BAM index..."

# samtools index creates a .bai file for random access to the BAM
# This is required by most downstream tools (featureCounts, IGV, etc.)
samtools index \
    -@ "${threads}" \
    "${output_bam}"

log "INFO" "BAM indexing completed"

# ==============================================================================
# OUTPUT VERIFICATION
# ==============================================================================
log "INFO" "Verifying outputs..."

[[ -s "${output_bam}" ]] || die "Output BAM file is empty or missing: ${output_bam}"
[[ -s "${output_bam}.bai" ]] || die "BAM index file is empty or missing: ${output_bam}.bai"

# Quick validation of BAM integrity
if ! samtools quickcheck "${output_bam}"; then
    die "BAM file failed integrity check: ${output_bam}"
fi

# ==============================================================================
# BAM STATISTICS
# ==============================================================================
log "INFO" "Generating BAM statistics..."

# Use samtools flagstat for alignment statistics
stats_file="${output_dir}/${sample_id}.flagstat.txt"
samtools flagstat \
    -@ "${threads}" \
    "${output_bam}" > "${stats_file}"

# Parse and log key statistics
if [[ -f "${stats_file}" ]]; then
    total_reads=$(grep "in total" "${stats_file}" | cut -d' ' -f1)
    mapped_reads=$(grep "mapped (" "${stats_file}" | head -1 | cut -d' ' -f1)
    paired_reads=$(grep "properly paired" "${stats_file}" | cut -d' ' -f1)

    log "INFO" "BAM statistics:"
    log "INFO" "  Total reads:      ${total_reads:-N/A}"
    log "INFO" "  Mapped reads:     ${mapped_reads:-N/A}"
    log "INFO" "  Properly paired:  ${paired_reads:-N/A}"
fi

# ==============================================================================
# COMPLETION
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "BAM processing completed for: ${sample_id}"
log "INFO" "Output BAM: ${output_bam}"
log "INFO" "=========================================="

exit 0
