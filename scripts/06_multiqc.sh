#!/usr/bin/env bash
#
# 06_multiqc.sh - Generate comprehensive QC report using MultiQC
#
# Description:
#   This script aggregates quality control metrics from multiple tools into a
#   single interactive HTML report using MultiQC. It searches the results
#   directory for outputs from FastP, STAR, RSeQC, Qualimap, featureCounts,
#   and other supported tools, presenting them in a unified format for easy
#   comparison across samples.
#
# Usage:
#   ./06_multiqc.sh -i <input_dir> -o <output_dir> [-t title]
#
# Arguments:
#   -i, --input         Path to directory containing QC outputs (required)
#   -o, --output-dir    Output directory for MultiQC report (required)
#   -t, --title         Report title (default: "RNA-seq QC Report")
#   -c, --config        Path to MultiQC config file (optional)
#   -f, --force         Overwrite existing reports
#   -h, --help          Display this help message
#
# Dependencies:
#   - MultiQC >= 1.32 (QC report aggregation)
#
# Output:
#   - multiqc_report.html         Interactive HTML report
#   - multiqc_data/               Directory containing parsed data
#
# Example:
#   ./06_multiqc.sh -i results/ -o results/06_multiqc/ -t "Project XYZ QC"
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
readonly DEFAULT_TITLE="RNA-seq QC Report"

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
Usage: ${SCRIPT_NAME} -i <input_dir> -o <output_dir> [OPTIONS]

Generate comprehensive QC report using MultiQC.

Required Arguments:
  -i, --input         Path to directory containing QC outputs (will search recursively)
  -o, --output-dir    Output directory for MultiQC report

Optional Arguments:
  -t, --title         Report title (default: "${DEFAULT_TITLE}")
  -c, --config        Path to MultiQC config file
  -f, --force         Overwrite existing reports
  -h, --help          Display this help message

Supported Input Formats:
  - FastP (fastp.json)
  - STAR (Log.final.out)
  - RSeQC (various outputs)
  - Qualimap (rnaseq_qc_results.txt)
  - featureCounts (*.summary)
  - samtools (flagstat, stats)

Examples:
  ${SCRIPT_NAME} -i results/ -o results/06_multiqc/
  ${SCRIPT_NAME} -i results/ -o results/06_multiqc/ -t "My Project" -f

EOF
    exit 0
}

check_command() {
    local cmd="$1"
    if ! command -v "${cmd}" &> /dev/null; then
        die "Required command '${cmd}' not found in PATH"
    fi
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

input_dir=""
output_dir=""
title="${DEFAULT_TITLE}"
config_file=""
force_flag=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input_dir="$2"
            shift 2
            ;;
        -o|--output-dir)
            output_dir="$2"
            shift 2
            ;;
        -t|--title)
            title="$2"
            shift 2
            ;;
        -c|--config)
            config_file="$2"
            shift 2
            ;;
        -f|--force)
            force_flag="--force"
            shift
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

[[ -z "${input_dir}" ]] && die "Input directory is required (-i)"
[[ -z "${output_dir}" ]] && die "Output directory is required (-o)"
[[ -d "${input_dir}" ]] || die "Input directory not found: ${input_dir}"

if [[ -n "${config_file}" ]] && [[ ! -f "${config_file}" ]]; then
    die "Config file not found: ${config_file}"
fi

check_command "multiqc"
mkdir -p "${output_dir}"

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Generating MultiQC Report"
log "INFO" "=========================================="

# Log MultiQC version
multiqc_version="$(multiqc --version 2>&1)"
log "INFO" "Using ${multiqc_version}"

log "INFO" "Scanning directory: ${input_dir}"
log "INFO" "Report title: ${title}"

# ==============================================================================
# RUN MULTIQC
# ==============================================================================
# MultiQC automatically detects and parses outputs from supported tools.
# The following tools from our pipeline are supported:
#   - FastP: Read filtering statistics
#   - STAR: Alignment statistics
#   - RSeQC: RNA-seq QC metrics
#   - Qualimap: General alignment QC
#   - featureCounts: Read counting summary
#   - samtools: flagstat and stats output

# Build MultiQC command
multiqc_cmd=(
    multiqc
    "${input_dir}"
    --outdir "${output_dir}"
    --title "${title}"
    # Create clean output structure
    --filename "multiqc_report"
    # Export data in multiple formats
    --export
    # Include plots data for downstream analysis
    --data-format json
    # Zip data directory to save space
    --zip-data-dir
    # Verbose output for debugging
    --verbose
)

# Add force flag if specified
if [[ -n "${force_flag}" ]]; then
    multiqc_cmd+=(--force)
fi

# Add config file if specified
if [[ -n "${config_file}" ]]; then
    multiqc_cmd+=(--config "${config_file}")
fi

log "INFO" "Running MultiQC..."

# Execute MultiQC
if "${multiqc_cmd[@]}"; then
    log "INFO" "MultiQC completed successfully"
else
    die "MultiQC failed"
fi

# ==============================================================================
# OUTPUT VERIFICATION
# ==============================================================================
log "INFO" "Verifying outputs..."

report_file="${output_dir}/multiqc_report.html"
[[ -f "${report_file}" ]] || die "MultiQC report not generated: ${report_file}"

# Log detected modules
log "INFO" "Detected data from the following tools:"
if [[ -f "${output_dir}/multiqc_data/multiqc_sources.txt" ]]; then
    # Extract unique module names from sources file
    cut -f1 "${output_dir}/multiqc_data/multiqc_sources.txt" | sort -u | while read -r module; do
        log "INFO" "  - ${module}"
    done
fi

# ==============================================================================
# COMPLETION
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "MultiQC report generated successfully"
log "INFO" "Report: ${report_file}"
log "INFO" "Data:   ${output_dir}/multiqc_data/"
log "INFO" "=========================================="

exit 0
