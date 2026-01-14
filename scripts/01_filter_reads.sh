#!/usr/bin/env bash
#
# 01_filter_reads.sh - Filter and trim raw sequencing reads using FastP
#
# Description:
#   This script performs quality control filtering on raw FASTQ files using FastP.
#   It applies poly-X tail trimming, 3' quality-based trimming, and filters reads
#   based on minimum quality score and length requirements. This is typically the
#   first processing step after demultiplexing.
#
# Usage:
#   ./01_filter_reads.sh -i <input_r1> -o <output_dir> [-p <input_r2>] [-t threads]
#
# Arguments:
#   -i, --input-r1      Path to input FASTQ file (R1, required)
#   -p, --input-r2      Path to input FASTQ file (R2, for paired-end)
#   -o, --output-dir    Output directory for filtered reads (required)
#   -s, --sample-id     Sample identifier (default: derived from input filename)
#   -t, --threads       Number of threads (default: 8)
#   -q, --min-quality   Minimum Phred quality score (default: 15)
#   -l, --min-length    Minimum read length after trimming (default: 50)
#   -h, --help          Display this help message
#
# Dependencies:
#   - fastp >= 0.24.0 (read filtering and QC)
#
# Output:
#   - {sample}_R1.filtered.fastq.gz  Filtered R1 reads
#   - {sample}_R2.filtered.fastq.gz  Filtered R2 reads (if paired-end)
#   - {sample}_fastp.html            HTML QC report
#   - {sample}_fastp.json            JSON QC metrics
#
# Example:
#   # Single-end reads
#   ./01_filter_reads.sh -i sample_R1.fastq.gz -o results/01_filtered/
#
#   # Paired-end reads
#   ./01_filter_reads.sh -i sample_R1.fastq.gz -p sample_R2.fastq.gz \
#       -o results/01_filtered/ -t 16
#
# Author: Plasmidsaurus RNA-seq Pipeline
# Date: 2024
# Version: 1.0.0

# ==============================================================================
# SHELL OPTIONS
# ==============================================================================
# Exit immediately on error (-e), treat unset variables as errors (-u),
# and ensure pipeline failures are caught (-o pipefail)
set -euo pipefail

# Enable debug mode if DEBUG environment variable is set
[[ "${DEBUG:-}" == "true" ]] && set -x

# ==============================================================================
# CONSTANTS
# ==============================================================================
readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Default parameter values
readonly DEFAULT_THREADS=8
readonly DEFAULT_MIN_QUALITY=15
readonly DEFAULT_MIN_LENGTH=50

# Tool version requirements
readonly FASTP_MIN_VERSION="0.24.0"

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#######################################
# Print timestamped log message to stderr.
#
# Arguments:
#   $1 - Log level (INFO, WARN, ERROR)
#   $2 - Log message
# Outputs:
#   Writes timestamped message to stderr
#######################################
log() {
    local level="$1"
    shift
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level}] $*" >&2
}

#######################################
# Print error message and exit with non-zero status.
#
# Arguments:
#   $1 - Error message
# Returns:
#   Exits with status 1
#######################################
die() {
    log "ERROR" "$*"
    exit 1
}

#######################################
# Display usage information and exit.
#
# Outputs:
#   Writes usage information to stdout
# Returns:
#   Exits with status 0
#######################################
usage() {
    cat << EOF
Usage: ${SCRIPT_NAME} -i <input_r1> -o <output_dir> [OPTIONS]

Filter and trim raw sequencing reads using FastP.

Required Arguments:
  -i, --input-r1      Path to input FASTQ file (R1)
  -o, --output-dir    Output directory for filtered reads

Optional Arguments:
  -p, --input-r2      Path to input FASTQ file (R2, for paired-end)
  -s, --sample-id     Sample identifier (default: derived from input filename)
  -t, --threads       Number of threads (default: ${DEFAULT_THREADS})
  -q, --min-quality   Minimum Phred quality score (default: ${DEFAULT_MIN_QUALITY})
  -l, --min-length    Minimum read length after trimming (default: ${DEFAULT_MIN_LENGTH})
  -h, --help          Display this help message

Examples:
  # Single-end reads
  ${SCRIPT_NAME} -i sample_R1.fastq.gz -o results/01_filtered/

  # Paired-end reads with custom parameters
  ${SCRIPT_NAME} -i sample_R1.fastq.gz -p sample_R2.fastq.gz \\
      -o results/01_filtered/ -t 16 -q 20 -l 75

EOF
    exit 0
}

#######################################
# Check if a required command is available.
#
# Arguments:
#   $1 - Command name to check
# Returns:
#   0 if command exists, exits with error otherwise
#######################################
check_command() {
    local cmd="$1"
    if ! command -v "${cmd}" &> /dev/null; then
        die "Required command '${cmd}' not found in PATH"
    fi
}

#######################################
# Extract sample ID from FASTQ filename.
#
# Removes common suffixes like _R1, _1, .fastq, .gz, etc.
#
# Arguments:
#   $1 - FASTQ filename
# Outputs:
#   Writes extracted sample ID to stdout
#######################################
extract_sample_id() {
    local filename="$1"
    # Remove directory path
    filename="$(basename "${filename}")"
    # Remove common suffixes
    filename="${filename%.gz}"
    filename="${filename%.fastq}"
    filename="${filename%.fq}"
    filename="${filename%_R1}"
    filename="${filename%_R2}"
    filename="${filename%_1}"
    filename="${filename%_2}"
    echo "${filename}"
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

# Initialize variables with defaults
input_r1=""
input_r2=""
output_dir=""
sample_id=""
threads="${DEFAULT_THREADS}"
min_quality="${DEFAULT_MIN_QUALITY}"
min_length="${DEFAULT_MIN_LENGTH}"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-r1)
            input_r1="$2"
            shift 2
            ;;
        -p|--input-r2)
            input_r2="$2"
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
        -q|--min-quality)
            min_quality="$2"
            shift 2
            ;;
        -l|--min-length)
            min_length="$2"
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

# Check required arguments
[[ -z "${input_r1}" ]] && die "Input R1 file is required (-i)"
[[ -z "${output_dir}" ]] && die "Output directory is required (-o)"

# Validate input files exist
[[ -f "${input_r1}" ]] || die "Input R1 file not found: ${input_r1}"
if [[ -n "${input_r2}" ]]; then
    [[ -f "${input_r2}" ]] || die "Input R2 file not found: ${input_r2}"
fi

# Derive sample ID if not provided
if [[ -z "${sample_id}" ]]; then
    sample_id="$(extract_sample_id "${input_r1}")"
    log "INFO" "Derived sample ID: ${sample_id}"
fi

# Check dependencies
check_command "fastp"

# Create output directory if it doesn't exist
mkdir -p "${output_dir}"

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Starting read filtering for: ${sample_id}"
log "INFO" "=========================================="

# Log tool version for reproducibility
fastp_version="$(fastp --version 2>&1 | head -1)"
log "INFO" "Using ${fastp_version}"

# Define output file paths
output_r1="${output_dir}/${sample_id}_R1.filtered.fastq.gz"
output_html="${output_dir}/${sample_id}_fastp.html"
output_json="${output_dir}/${sample_id}_fastp.json"

# Build the fastp command
# Using an array ensures proper handling of arguments with spaces
fastp_cmd=(
    fastp
    --in1 "${input_r1}"
    --out1 "${output_r1}"
    --html "${output_html}"
    --json "${output_json}"
    --thread "${threads}"
    # Quality filtering: minimum Phred score for qualified bases
    --qualified_quality_phred "${min_quality}"
    # Length filtering: discard reads shorter than this after trimming
    --length_required "${min_length}"
    # Poly-X trimming: remove poly-X tails (common in RNA-seq)
    # This addresses issues with poly-A tails and other homopolymer artifacts
    --trim_poly_x
    # 3' quality trimming: sliding window approach
    # Window size of 4 bases, trim if mean quality drops below threshold
    --cut_tail
    --cut_tail_window_size 4
    --cut_tail_mean_quality "${min_quality}"
    # Disable adapter auto-detection logging noise
    --dont_eval_duplication
)

# Add paired-end specific options if R2 is provided
if [[ -n "${input_r2}" ]]; then
    output_r2="${output_dir}/${sample_id}_R2.filtered.fastq.gz"
    fastp_cmd+=(
        --in2 "${input_r2}"
        --out2 "${output_r2}"
        # Detect and trim adapters for paired-end reads
        --detect_adapter_for_pe
    )
    log "INFO" "Running in paired-end mode"
else
    log "INFO" "Running in single-end mode"
fi

# Log the command being executed
log "INFO" "Executing fastp command..."
log "INFO" "Parameters: min_quality=${min_quality}, min_length=${min_length}, threads=${threads}"

# Execute fastp
if "${fastp_cmd[@]}"; then
    log "INFO" "FastP filtering completed successfully"
else
    die "FastP filtering failed for sample: ${sample_id}"
fi

# ==============================================================================
# OUTPUT VERIFICATION
# ==============================================================================
log "INFO" "Verifying outputs..."

# Check that output files were created and are not empty
[[ -s "${output_r1}" ]] || die "Output R1 file is empty or missing: ${output_r1}"
[[ -s "${output_json}" ]] || die "JSON report is empty or missing: ${output_json}"

if [[ -n "${input_r2}" ]]; then
    [[ -s "${output_r2}" ]] || die "Output R2 file is empty or missing: ${output_r2}"
fi

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
# Parse key metrics from the JSON report for logging
if command -v jq &> /dev/null; then
    log "INFO" "Filtering summary:"

    total_reads_before=$(jq '.summary.before_filtering.total_reads' "${output_json}")
    total_reads_after=$(jq '.summary.after_filtering.total_reads' "${output_json}")
    q30_before=$(jq '.summary.before_filtering.q30_rate' "${output_json}")
    q30_after=$(jq '.summary.after_filtering.q30_rate' "${output_json}")

    # Calculate retention rate
    if [[ "${total_reads_before}" -gt 0 ]]; then
        retention_rate=$(echo "scale=2; ${total_reads_after} * 100 / ${total_reads_before}" | bc)
        log "INFO" "  Total reads before: ${total_reads_before}"
        log "INFO" "  Total reads after:  ${total_reads_after}"
        log "INFO" "  Retention rate:     ${retention_rate}%"
        log "INFO" "  Q30 before:         ${q30_before}"
        log "INFO" "  Q30 after:          ${q30_after}"
    fi
else
    log "INFO" "Install 'jq' for detailed filtering statistics"
fi

# ==============================================================================
# COMPLETION
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Read filtering completed for: ${sample_id}"
log "INFO" "Output directory: ${output_dir}"
log "INFO" "=========================================="

exit 0
