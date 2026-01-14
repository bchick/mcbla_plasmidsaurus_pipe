#!/usr/bin/env bash
#
# run_pipeline.sh - Main orchestrator for the Plasmidsaurus RNA-seq pipeline
#
# Description:
#   This script orchestrates the complete RNA-seq analysis pipeline from raw
#   FASTQ files or deduplicated BAM files through differential expression and
#   functional enrichment. It manages sample processing, handles dependencies
#   between steps, and provides options for resuming from specific stages.
#
# Usage:
#   ./run_pipeline.sh -i <input_dir> -o <output_dir> [-g genome | -c config] [OPTIONS]
#
# Arguments:
#   -i, --input         Directory containing input files (required)
#   -o, --output        Output directory for results (required)
#   -g, --genome        Reference genome: 'hg38' or 'mm10' (auto-selects config)
#   -c, --config        Path to YAML configuration file (alternative to -g)
#   -y, --type          Input data type: 'fastq' or 'bam' (default: fastq)
#   -s, --start-step    Step to start from (1-10, default: 1)
#   -e, --end-step      Step to end at (1-10, default: 10)
#   -t, --threads       Number of threads (default: 8)
#   -m, --metadata      Path to sample metadata TSV file (required for steps 8+)
#   -d, --dry-run       Print commands without executing
#   -h, --help          Display this help message
#
# Pipeline Steps:
#   1  - Read filtering (FastP)
#   2  - Alignment (STAR)
#   3  - BAM sorting (samtools)
#   4  - Deduplication (UMICollapse)
#   5  - Mapping QC (RSeQC/Qualimap)    <- BAM input starts here
#   6  - QC Report (MultiQC)
#   7  - Quantification (featureCounts)
#   8  - Normalization (edgeR TMM)
#   9  - Differential Expression (edgeR)
#   10 - Functional Enrichment (GSEApy)
#
# Dependencies:
#   All tools listed in README.md must be installed and in PATH
#
# Example:
#   # Run full pipeline from FASTQ (human)
#   ./run_pipeline.sh -i data/fastq/ -o results/ -g hg38 -m config/samples.tsv
#
#   # Run from deduplicated BAM files (mouse)
#   ./run_pipeline.sh -i data/bam/ -o results/ -y bam -g mm10 -m config/samples.tsv
#
#   # Run only QC steps
#   ./run_pipeline.sh -i data/fastq/ -o results/ -g hg38 -s 1 -e 6
#
# Author: Plasmidsaurus RNA-seq Pipeline
# Date: 2024
# Version: 1.1.0

# ==============================================================================
# SHELL OPTIONS
# ==============================================================================
set -euo pipefail
[[ "${DEBUG:-}" == "true" ]] && set -x

# ==============================================================================
# CONSTANTS
# ==============================================================================
readonly SCRIPT_NAME="$(basename "$0")"
readonly SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
readonly SCRIPTS_DIR="${SCRIPT_DIR}/scripts"

# Pipeline step names for logging
readonly STEP_NAMES=(
    ""  # Index 0 unused
    "Read Filtering (FastP)"
    "Alignment (STAR)"
    "BAM Sorting (samtools)"
    "Deduplication (UMICollapse)"
    "Mapping QC (RSeQC/Qualimap)"
    "QC Report (MultiQC)"
    "Quantification (featureCounts)"
    "Normalization (edgeR)"
    "Differential Expression (edgeR)"
    "Functional Enrichment (GSEApy)"
)

# Default values
readonly DEFAULT_THREADS=8
readonly DEFAULT_START_STEP=1
readonly DEFAULT_END_STEP=10

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
Usage: ${SCRIPT_NAME} -i <input_dir> -o <output_dir> [-g genome | -c config] [OPTIONS]

Plasmidsaurus RNA-seq Analysis Pipeline

Required Arguments:
  -i, --input         Directory containing input files (FASTQ or BAM)
  -o, --output        Output directory for results
  -g, --genome        Reference genome: 'hg38' (human) or 'mm10' (mouse)
                      Auto-selects config file (alternative to -c)
  -c, --config        Path to YAML configuration file (alternative to -g)

Optional Arguments:
  -y, --type          Input data type: 'fastq' or 'bam' (default: fastq)
                      BAM input automatically starts at step 5 (Mapping QC)
  -s, --start-step    Step to start from (1-10, default: ${DEFAULT_START_STEP})
  -e, --end-step      Step to end at (1-10, default: ${DEFAULT_END_STEP})
  -t, --threads       Number of threads (default: ${DEFAULT_THREADS})
  -m, --metadata      Path to sample metadata TSV file (required for steps 8+)
  -d, --dry-run       Print commands without executing
  -h, --help          Display this help message

Pipeline Steps:
  1  - Read Filtering (FastP)
  2  - Alignment (STAR)
  3  - BAM Sorting (samtools)
  4  - Deduplication (UMICollapse)
  5  - Mapping QC (RSeQC/Qualimap)    <- BAM input starts here
  6  - QC Report (MultiQC)
  7  - Quantification (featureCounts)
  8  - Normalization (edgeR TMM)
  9  - Differential Expression (edgeR)
  10 - Functional Enrichment (GSEApy)

Examples:
  # Run full pipeline from FASTQ (human)
  ${SCRIPT_NAME} -i data/fastq/ -o results/ -g hg38 -m config/samples.tsv

  # Run from deduplicated BAM files (mouse)
  ${SCRIPT_NAME} -i data/bam/ -o results/ -y bam -g mm10 -m config/samples.tsv

  # Run only preprocessing (steps 1-4)
  ${SCRIPT_NAME} -i data/fastq/ -o results/ -g hg38 -s 1 -e 4

  # Use custom config file
  ${SCRIPT_NAME} -i data/fastq/ -o results/ -c config/config.yaml -m config/samples.tsv

EOF
    exit 0
}

#######################################
# Parse YAML config file (basic parser for simple YAML)
#
# Arguments:
#   $1 - Config file path
#   $2 - Key to retrieve (dot notation: section.key)
# Outputs:
#   Value from config file
#######################################
get_config() {
    local config_file="$1"
    local key="$2"

    # Simple YAML parser using grep/sed
    # For production use, consider using yq or a proper YAML parser
    local value
    value=$(grep -E "^\s*${key##*.}:" "${config_file}" | head -1 | sed 's/.*:\s*//' | sed 's/\s*#.*//' | tr -d '"' | tr -d "'")

    echo "${value}"
}

#######################################
# Print step header for visual separation in logs
#
# Arguments:
#   $1 - Step number
#   $2 - Step name
#######################################
print_step_header() {
    local step_num="$1"
    local step_name="$2"

    log "INFO" ""
    log "INFO" "============================================================"
    log "INFO" "STEP ${step_num}: ${step_name}"
    log "INFO" "============================================================"
    log "INFO" ""
}

#######################################
# Find FASTQ files in a directory
#
# Arguments:
#   $1 - Directory to search
# Outputs:
#   List of FASTQ files (R1 only for paired-end detection)
#######################################
find_fastq_files() {
    local input_dir="$1"

    # Find R1 files (will pair with R2 later)
    find "${input_dir}" -type f \( -name "*_R1*.fastq.gz" -o -name "*_R1*.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort
}

#######################################
# Get R2 file path from R1 file path
#
# Arguments:
#   $1 - R1 file path
# Outputs:
#   R2 file path (or empty if not found/single-end)
#######################################
get_r2_file() {
    local r1_file="$1"
    local r2_file

    # Try common R2 naming patterns
    # Only return if substitution actually changed the path AND file exists
    r2_file="${r1_file/_R1/_R2}"
    [[ "${r2_file}" != "${r1_file}" ]] && [[ -f "${r2_file}" ]] && echo "${r2_file}" && return

    r2_file="${r1_file/_1.fastq/_2.fastq}"
    [[ "${r2_file}" != "${r1_file}" ]] && [[ -f "${r2_file}" ]] && echo "${r2_file}" && return

    r2_file="${r1_file/_1.fq/_2.fq}"
    [[ "${r2_file}" != "${r1_file}" ]] && [[ -f "${r2_file}" ]] && echo "${r2_file}" && return

    # No R2 found (single-end)
    echo ""
}

#######################################
# Extract sample ID from FASTQ filename
#
# Arguments:
#   $1 - FASTQ file path
# Outputs:
#   Sample ID
#######################################
extract_sample_id() {
    local filepath="$1"
    local filename

    filename="$(basename "${filepath}")"
    filename="${filename%.gz}"
    filename="${filename%.fastq}"
    filename="${filename%.fq}"
    filename="${filename%_R1}"
    filename="${filename%_R2}"
    filename="${filename%_1}"
    filename="${filename%_2}"

    echo "${filename}"
}

#######################################
# Find BAM files in a directory
#
# Arguments:
#   $1 - Directory to search
# Outputs:
#   List of BAM file paths (excluding .bai index files)
#######################################
find_bam_files() {
    local input_dir="$1"

    # Find all BAM files, excluding index files
    find "${input_dir}" -type f -name "*.bam" ! -name "*.bai" | sort
}

#######################################
# Extract sample ID from BAM filename
#
# Strips common suffixes from BAM filenames to extract the sample ID.
# Handles various naming conventions from different pipeline outputs.
#
# Arguments:
#   $1 - BAM file path
# Outputs:
#   Sample ID
#######################################
extract_sample_id_from_bam() {
    local filepath="$1"
    local filename

    filename="$(basename "${filepath}" .bam)"
    # Strip common pipeline output suffixes
    filename="${filename%_dedup-mapped-reads}"
    filename="${filename%.dedup}"
    filename="${filename%.sorted}"
    filename="${filename%_Aligned.sortedByCoord.out}"
    filename="${filename%_Aligned.out}"

    echo "${filename}"
}

#######################################
# Find all FASTQ files (single-end) in a directory
#
# Unlike find_fastq_files(), this finds ALL fastq.gz files regardless
# of R1/R2 naming conventions - useful for single-end data.
#
# Arguments:
#   $1 - Directory to search
# Outputs:
#   List of all FASTQ file paths
#######################################
find_all_fastq_files() {
    local input_dir="$1"

    find "${input_dir}" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | sort
}

#######################################
# Extract sample ID from single-end FASTQ filename
#
# Simply strips the .fastq.gz extension without removing R1/R2 suffixes.
#
# Arguments:
#   $1 - FASTQ file path
# Outputs:
#   Sample ID (filename without extension)
#######################################
extract_sample_id_single_end() {
    local filepath="$1"
    local filename

    filename="$(basename "${filepath}")"
    filename="${filename%.gz}"
    filename="${filename%.fastq}"
    filename="${filename%.fq}"

    echo "${filename}"
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

input_dir=""
output_dir=""
config_file=""
metadata_file=""
start_step="${DEFAULT_START_STEP}"
end_step="${DEFAULT_END_STEP}"
threads="${DEFAULT_THREADS}"
dry_run=false
input_type="fastq"    # Input data type: 'fastq' or 'bam'
genome=""             # Genome shorthand: 'hg38' or 'mm10'

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input_dir="$2"
            shift 2
            ;;
        -o|--output)
            output_dir="$2"
            shift 2
            ;;
        -c|--config)
            config_file="$2"
            shift 2
            ;;
        -m|--metadata)
            metadata_file="$2"
            shift 2
            ;;
        -s|--start-step)
            start_step="$2"
            shift 2
            ;;
        -e|--end-step)
            end_step="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -d|--dry-run)
            dry_run=true
            shift
            ;;
        -y|--type)
            input_type="$2"
            shift 2
            ;;
        -g|--genome)
            genome="$2"
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

log "INFO" "============================================================"
log "INFO" "PLASMIDSAURUS RNA-SEQ ANALYSIS PIPELINE"
log "INFO" "============================================================"
log "INFO" ""

# Validate input type
[[ "${input_type}" == "fastq" ]] || [[ "${input_type}" == "bam" ]] || \
    die "Invalid input type: ${input_type}. Use 'fastq' or 'bam'"

# Auto-select config file based on genome if provided
if [[ -n "${genome}" ]]; then
    case "${genome}" in
        hg38) config_file="${SCRIPT_DIR}/config/config.human.yaml" ;;
        mm10) config_file="${SCRIPT_DIR}/config/config.mouse.yaml" ;;
        *) die "Unknown genome: ${genome}. Use 'hg38' or 'mm10'" ;;
    esac
    log "INFO" "Auto-selected config for genome ${genome}: ${config_file}"
fi

# Check required arguments
[[ -z "${input_dir}" ]] && die "Input directory is required (-i)"
[[ -z "${output_dir}" ]] && die "Output directory is required (-o)"
[[ -z "${config_file}" ]] && die "Configuration file is required (-c or -g/--genome)"

# Validate paths
[[ -d "${input_dir}" ]] || die "Input directory not found: ${input_dir}"
[[ -f "${config_file}" ]] || die "Config file not found: ${config_file}"

# Validate step range
[[ "${start_step}" -ge 1 ]] && [[ "${start_step}" -le 10 ]] || die "Start step must be 1-10"
[[ "${end_step}" -ge 1 ]] && [[ "${end_step}" -le 10 ]] || die "End step must be 1-10"
[[ "${start_step}" -le "${end_step}" ]] || die "Start step must be <= end step"

# Check metadata file for downstream steps
if [[ "${end_step}" -ge 8 ]] && [[ -z "${metadata_file}" ]]; then
    die "Metadata file is required for steps 8+ (-m)"
fi
if [[ -n "${metadata_file}" ]] && [[ ! -f "${metadata_file}" ]]; then
    die "Metadata file not found: ${metadata_file}"
fi

# Create output directories
mkdir -p "${output_dir}"/{01_filtered,02_aligned,03_sorted,04_dedup,05_qc,06_multiqc,07_counts,08_normalization,09_de,10_enrichment}
mkdir -p "${output_dir}/logs"

# Log configuration
log "INFO" "Configuration:"
log "INFO" "  Input directory:   ${input_dir}"
log "INFO" "  Input type:        ${input_type}"
log "INFO" "  Output directory:  ${output_dir}"
log "INFO" "  Genome:            ${genome:-N/A}"
log "INFO" "  Config file:       ${config_file}"
log "INFO" "  Metadata file:     ${metadata_file:-N/A}"
log "INFO" "  Steps to run:      ${start_step} - ${end_step}"
log "INFO" "  Threads:           ${threads}"
log "INFO" "  Dry run:           ${dry_run}"
log "INFO" ""

# ==============================================================================
# LOAD CONFIGURATION
# ==============================================================================

log "INFO" "Loading configuration..."

# Extract key paths from config
genome_fasta=$(get_config "${config_file}" "genome_fasta")
gtf_file=$(get_config "${config_file}" "gtf")
star_index=$(get_config "${config_file}" "star_index")
gene_model_bed=$(get_config "${config_file}" "gene_model_bed")
organism=$(get_config "${config_file}" "organism")

# Use defaults if not specified
organism="${organism:-human}"

log "INFO" "Reference files:"
log "INFO" "  Genome FASTA:  ${genome_fasta:-Not specified}"
log "INFO" "  GTF file:      ${gtf_file:-Not specified}"
log "INFO" "  STAR index:    ${star_index:-Not specified}"
log "INFO" "  Gene BED:      ${gene_model_bed:-Not specified}"
log "INFO" "  Organism:      ${organism}"

# ==============================================================================
# DISCOVER SAMPLES
# ==============================================================================

log "INFO" "Discovering samples..."
log "INFO" "  Input type: ${input_type}"

# Arrays to store input files and sample IDs
declare -a input_files
declare -a sample_ids

if [[ "${input_type}" == "fastq" ]]; then
    # Try paired-end naming first, fall back to all FASTQ files (single-end)
    mapfile -t input_files < <(find_fastq_files "${input_dir}")

    if [[ ${#input_files[@]} -eq 0 ]]; then
        # No paired-end files found, try single-end (all FASTQ files)
        log "INFO" "No paired-end FASTQ files found, searching for single-end..."
        mapfile -t input_files < <(find_all_fastq_files "${input_dir}")

        if [[ ${#input_files[@]} -eq 0 ]]; then
            die "No FASTQ files found in ${input_dir}"
        fi

        # Extract sample IDs for single-end data
        for fq in "${input_files[@]}"; do
            sample_ids+=("$(extract_sample_id_single_end "${fq}")")
        done
    else
        # Extract sample IDs for paired-end data
        for fq in "${input_files[@]}"; do
            sample_ids+=("$(extract_sample_id "${fq}")")
        done
    fi

    # For backward compatibility, also populate fastq_files array
    fastq_files=("${input_files[@]}")

elif [[ "${input_type}" == "bam" ]]; then
    mapfile -t input_files < <(find_bam_files "${input_dir}")

    if [[ ${#input_files[@]} -eq 0 ]]; then
        die "No BAM files found in ${input_dir}"
    fi

    # Extract sample IDs from BAM filenames
    for bam in "${input_files[@]}"; do
        sample_ids+=("$(extract_sample_id_from_bam "${bam}")")
    done

    # BAM input should start at step 5 (mapping QC) at earliest
    if [[ "${start_step}" -lt 5 ]]; then
        log "INFO" "BAM input detected - adjusting start step from ${start_step} to 5 (Mapping QC)"
        start_step=5
    fi
fi

log "INFO" "Found ${#sample_ids[@]} samples:"
for i in "${!sample_ids[@]}"; do
    log "INFO" "  - ${sample_ids[i]} (${input_files[i]})"
done

# ==============================================================================
# PIPELINE EXECUTION
# ==============================================================================

# Helper function to run or print command
run_cmd() {
    if [[ "${dry_run}" == true ]]; then
        echo "[DRY-RUN] $*"
    else
        "$@"
    fi
}

# Track pipeline start time
pipeline_start=$(date +%s)

# ------------------------------------------------------------------------------
# STEP 1: Read Filtering (FastP)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 1 ]] && [[ "${end_step}" -ge 1 ]]; then
    print_step_header 1 "${STEP_NAMES[1]}"

    for r1_file in "${fastq_files[@]}"; do
        sample_id=$(extract_sample_id "${r1_file}")
        r2_file=$(get_r2_file "${r1_file}")

        log "INFO" "Processing sample: ${sample_id}"

        cmd=("${SCRIPTS_DIR}/01_filter_reads.sh"
            -i "${r1_file}"
            -o "${output_dir}/01_filtered"
            -s "${sample_id}"
            -t "${threads}")

        [[ -n "${r2_file}" ]] && cmd+=(-p "${r2_file}")

        run_cmd "${cmd[@]}" 2>&1 | tee -a "${output_dir}/logs/01_filter_${sample_id}.log"
    done

    log "INFO" "Step 1 completed"
fi

# ------------------------------------------------------------------------------
# STEP 2: Alignment (STAR)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 2 ]] && [[ "${end_step}" -ge 2 ]]; then
    print_step_header 2 "${STEP_NAMES[2]}"

    [[ -z "${star_index}" ]] && die "STAR index not specified in config"
    [[ -d "${star_index}" ]] || die "STAR index directory not found: ${star_index}"

    for sample_id in "${sample_ids[@]}"; do
        log "INFO" "Processing sample: ${sample_id}"

        r1_filtered="${output_dir}/01_filtered/${sample_id}_R1.filtered.fastq.gz"
        r2_filtered="${output_dir}/01_filtered/${sample_id}_R2.filtered.fastq.gz"

        cmd=("${SCRIPTS_DIR}/02_align.sh"
            -i "${r1_filtered}"
            -g "${star_index}"
            -o "${output_dir}/02_aligned"
            -s "${sample_id}"
            -t "${threads}")

        [[ -f "${r2_filtered}" ]] && cmd+=(-p "${r2_filtered}")

        run_cmd "${cmd[@]}" 2>&1 | tee -a "${output_dir}/logs/02_align_${sample_id}.log"
    done

    log "INFO" "Step 2 completed"
fi

# ------------------------------------------------------------------------------
# STEP 3: BAM Sorting (samtools)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 3 ]] && [[ "${end_step}" -ge 3 ]]; then
    print_step_header 3 "${STEP_NAMES[3]}"

    for sample_id in "${sample_ids[@]}"; do
        log "INFO" "Processing sample: ${sample_id}"

        input_bam="${output_dir}/02_aligned/${sample_id}_Aligned.sortedByCoord.out.bam"

        run_cmd "${SCRIPTS_DIR}/03_sort_bam.sh" \
            -i "${input_bam}" \
            -o "${output_dir}/03_sorted" \
            -s "${sample_id}" \
            -t "${threads}" \
            2>&1 | tee -a "${output_dir}/logs/03_sort_${sample_id}.log"
    done

    log "INFO" "Step 3 completed"
fi

# ------------------------------------------------------------------------------
# STEP 4: Deduplication (UMICollapse)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 4 ]] && [[ "${end_step}" -ge 4 ]]; then
    print_step_header 4 "${STEP_NAMES[4]}"

    for sample_id in "${sample_ids[@]}"; do
        log "INFO" "Processing sample: ${sample_id}"

        input_bam="${output_dir}/03_sorted/${sample_id}.sorted.bam"

        run_cmd "${SCRIPTS_DIR}/04_dedup.sh" \
            -i "${input_bam}" \
            -o "${output_dir}/04_dedup" \
            -s "${sample_id}" \
            2>&1 | tee -a "${output_dir}/logs/04_dedup_${sample_id}.log"
    done

    log "INFO" "Step 4 completed"
fi

# ------------------------------------------------------------------------------
# STEP 5: Mapping QC (RSeQC/Qualimap)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 5 ]] && [[ "${end_step}" -ge 5 ]]; then
    print_step_header 5 "${STEP_NAMES[5]}"

    for i in "${!sample_ids[@]}"; do
        sample_id="${sample_ids[i]}"
        log "INFO" "Processing sample: ${sample_id}"

        # Use original input BAM if input_type is bam, otherwise use pipeline output
        if [[ "${input_type}" == "bam" ]]; then
            input_bam="${input_files[i]}"
        else
            input_bam="${output_dir}/04_dedup/${sample_id}.dedup.bam"
        fi

        cmd=("${SCRIPTS_DIR}/05_mapping_qc.sh"
            -i "${input_bam}"
            -o "${output_dir}/05_qc"
            -s "${sample_id}"
            -t "${threads}")

        [[ -n "${gene_model_bed}" ]] && cmd+=(-b "${gene_model_bed}")
        [[ -n "${gtf_file}" ]] && cmd+=(-g "${gtf_file}")

        run_cmd "${cmd[@]}" 2>&1 | tee -a "${output_dir}/logs/05_qc_${sample_id}.log"
    done

    log "INFO" "Step 5 completed"
fi

# ------------------------------------------------------------------------------
# STEP 6: QC Report (MultiQC)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 6 ]] && [[ "${end_step}" -ge 6 ]]; then
    print_step_header 6 "${STEP_NAMES[6]}"

    run_cmd "${SCRIPTS_DIR}/06_multiqc.sh" \
        -i "${output_dir}" \
        -o "${output_dir}/06_multiqc" \
        -t "Plasmidsaurus RNA-seq QC Report" \
        -f \
        2>&1 | tee -a "${output_dir}/logs/06_multiqc.log"

    log "INFO" "Step 6 completed"
fi

# ------------------------------------------------------------------------------
# STEP 7: Quantification (featureCounts)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 7 ]] && [[ "${end_step}" -ge 7 ]]; then
    print_step_header 7 "${STEP_NAMES[7]}"

    [[ -z "${gtf_file}" ]] && die "GTF file not specified in config"
    [[ -f "${gtf_file}" ]] || die "GTF file not found: ${gtf_file}"

    # Use original input directory for BAM input, otherwise use pipeline dedup output
    if [[ "${input_type}" == "bam" ]]; then
        bam_input_dir="${input_dir}"
    else
        bam_input_dir="${output_dir}/04_dedup"
    fi

    run_cmd "${SCRIPTS_DIR}/07_quantify.sh" \
        -i "${bam_input_dir}" \
        -g "${gtf_file}" \
        -o "${output_dir}/07_counts" \
        -s 2 \
        -t "${threads}" \
        2>&1 | tee -a "${output_dir}/logs/07_quantify.log"

    log "INFO" "Step 7 completed"
fi

# ------------------------------------------------------------------------------
# STEP 8: Normalization (edgeR TMM)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 8 ]] && [[ "${end_step}" -ge 8 ]]; then
    print_step_header 8 "${STEP_NAMES[8]}"

    run_cmd Rscript "${SCRIPTS_DIR}/08_normalize.R" \
        --counts "${output_dir}/07_counts/gene_counts.txt" \
        --metadata "${metadata_file}" \
        --output "${output_dir}/08_normalization" \
        2>&1 | tee -a "${output_dir}/logs/08_normalize.log"

    log "INFO" "Step 8 completed"
fi

# ------------------------------------------------------------------------------
# STEP 9: Differential Expression (edgeR)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 9 ]] && [[ "${end_step}" -ge 9 ]]; then
    print_step_header 9 "${STEP_NAMES[9]}"

    # Get contrast from config (default to a placeholder)
    contrast=$(get_config "${config_file}" "contrasts")
    contrast="${contrast:-treatment-control}"

    run_cmd Rscript "${SCRIPTS_DIR}/09_differential_expression.R" \
        --dge "${output_dir}/08_normalization/dge_normalized.rds" \
        --metadata "${metadata_file}" \
        --output "${output_dir}/09_de" \
        --contrast "${contrast}" \
        2>&1 | tee -a "${output_dir}/logs/09_de.log"

    log "INFO" "Step 9 completed"
fi

# ------------------------------------------------------------------------------
# STEP 10: Functional Enrichment (GSEApy)
# ------------------------------------------------------------------------------
if [[ "${start_step}" -le 10 ]] && [[ "${end_step}" -ge 10 ]]; then
    print_step_header 10 "${STEP_NAMES[10]}"

    # Find ranked gene files from DE analysis
    for rnk_file in "${output_dir}/09_de/"ranked_genes_*.rnk; do
        [[ -f "${rnk_file}" ]] || continue

        contrast_name=$(basename "${rnk_file}" .rnk | sed 's/ranked_genes_//')
        log "INFO" "Running enrichment for contrast: ${contrast_name}"

        run_cmd python3 "${SCRIPTS_DIR}/10_enrichment.py" \
            --ranked "${rnk_file}" \
            --output "${output_dir}/10_enrichment/${contrast_name}" \
            --organism "${organism}" \
            --gene-sets "MSigDB_Hallmark_2020" \
            2>&1 | tee -a "${output_dir}/logs/10_enrichment_${contrast_name}.log"
    done

    log "INFO" "Step 10 completed"
fi

# ==============================================================================
# PIPELINE COMPLETION
# ==============================================================================

pipeline_end=$(date +%s)
runtime=$((pipeline_end - pipeline_start))
runtime_min=$((runtime / 60))
runtime_sec=$((runtime % 60))

log "INFO" ""
log "INFO" "============================================================"
log "INFO" "PIPELINE COMPLETED SUCCESSFULLY"
log "INFO" "============================================================"
log "INFO" ""
log "INFO" "Total runtime: ${runtime_min} minutes ${runtime_sec} seconds"
log "INFO" "Output directory: ${output_dir}"
log "INFO" ""
log "INFO" "Key outputs:"
log "INFO" "  - QC Report:      ${output_dir}/06_multiqc/multiqc_report.html"
log "INFO" "  - Count Matrix:   ${output_dir}/07_counts/gene_counts.txt"
log "INFO" "  - DE Results:     ${output_dir}/09_de/"
log "INFO" "  - Enrichment:     ${output_dir}/10_enrichment/"
log "INFO" ""

exit 0
