#!/usr/bin/env bash
#
# 02_align.sh - Align filtered reads to reference genome using STAR
#
# Description:
#   This script performs splice-aware alignment of RNA-seq reads to a reference
#   genome using the STAR aligner. STAR is optimized for RNA-seq data and can
#   identify both canonical and novel splice junctions. The script applies
#   non-canonical splice junction filtering and outputs unmapped reads for
#   downstream analysis or troubleshooting.
#
# Usage:
#   ./02_align.sh -i <input_r1> -g <genome_index> -o <output_dir> [-p <input_r2>] [-t threads]
#
# Arguments:
#   -i, --input-r1      Path to filtered FASTQ file (R1, required)
#   -p, --input-r2      Path to filtered FASTQ file (R2, for paired-end)
#   -g, --genome-index  Path to STAR genome index directory (required)
#   -o, --output-dir    Output directory for alignments (required)
#   -s, --sample-id     Sample identifier (default: derived from input filename)
#   -t, --threads       Number of threads (default: 16)
#   -h, --help          Display this help message
#
# Dependencies:
#   - STAR >= 2.7.11 (splice-aware aligner)
#
# Output:
#   - {sample}_Aligned.sortedByCoord.out.bam  Coordinate-sorted BAM file
#   - {sample}_Log.final.out                  Alignment summary statistics
#   - {sample}_Log.out                        Detailed run log
#   - {sample}_SJ.out.tab                     Detected splice junctions
#   - {sample}_Unmapped.out.mate1             Unmapped R1 reads (FASTQ)
#   - {sample}_Unmapped.out.mate2             Unmapped R2 reads (FASTQ, if PE)
#
# Example:
#   ./02_align.sh -i sample_R1.filtered.fastq.gz -p sample_R2.filtered.fastq.gz \
#       -g /path/to/star_index -o results/02_aligned/ -t 32
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
readonly SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Default parameter values
readonly DEFAULT_THREADS=16

# STAR alignment parameters (optimized for standard RNA-seq)
# These values are based on ENCODE RNA-seq guidelines and Plasmidsaurus recommendations
readonly DEFAULT_MULTIMAP_MAX=20           # Max multi-mapping locations
readonly DEFAULT_MISMATCH_MAX=999          # Max mismatches (effectively unlimited)
readonly DEFAULT_MISMATCH_RATIO=0.04       # Max mismatch ratio to read length
readonly DEFAULT_SJ_OVERHANG_UNANNOTATED=8 # Min overhang for unannotated junctions
readonly DEFAULT_SJ_OVERHANG_ANNOTATED=1   # Min overhang for annotated junctions
readonly DEFAULT_INTRON_MIN=20             # Minimum intron size
readonly DEFAULT_INTRON_MAX=1000000        # Maximum intron size (1 Mb)

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
Usage: ${SCRIPT_NAME} -i <input_r1> -g <genome_index> -o <output_dir> [OPTIONS]

Align filtered reads to reference genome using STAR.

Required Arguments:
  -i, --input-r1      Path to filtered FASTQ file (R1)
  -g, --genome-index  Path to STAR genome index directory
  -o, --output-dir    Output directory for alignments

Optional Arguments:
  -p, --input-r2      Path to filtered FASTQ file (R2, for paired-end)
  -s, --sample-id     Sample identifier (default: derived from input filename)
  -t, --threads       Number of threads (default: ${DEFAULT_THREADS})
  -h, --help          Display this help message

Examples:
  # Paired-end alignment
  ${SCRIPT_NAME} -i sample_R1.filtered.fastq.gz -p sample_R2.filtered.fastq.gz \\
      -g /path/to/star_index -o results/02_aligned/ -t 32

  # Single-end alignment
  ${SCRIPT_NAME} -i sample_R1.filtered.fastq.gz \\
      -g /path/to/star_index -o results/02_aligned/

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
    filename="${filename%.gz}"
    filename="${filename%.fastq}"
    filename="${filename%.fq}"
    filename="${filename%.filtered}"
    filename="${filename%_R1}"
    filename="${filename%_R2}"
    filename="${filename%_1}"
    filename="${filename%_2}"
    echo "${filename}"
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

input_r1=""
input_r2=""
genome_index=""
output_dir=""
sample_id=""
threads="${DEFAULT_THREADS}"

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
        -g|--genome-index)
            genome_index="$2"
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

# Check required arguments
[[ -z "${input_r1}" ]] && die "Input R1 file is required (-i)"
[[ -z "${genome_index}" ]] && die "Genome index directory is required (-g)"
[[ -z "${output_dir}" ]] && die "Output directory is required (-o)"

# Validate input files exist
[[ -f "${input_r1}" ]] || die "Input R1 file not found: ${input_r1}"
if [[ -n "${input_r2}" ]]; then
    [[ -f "${input_r2}" ]] || die "Input R2 file not found: ${input_r2}"
fi

# Validate genome index directory exists and contains required files
[[ -d "${genome_index}" ]] || die "Genome index directory not found: ${genome_index}"
[[ -f "${genome_index}/SA" ]] || die "Invalid STAR index: SA file not found in ${genome_index}"

# Derive sample ID if not provided
if [[ -z "${sample_id}" ]]; then
    sample_id="$(extract_sample_id "${input_r1}")"
    log "INFO" "Derived sample ID: ${sample_id}"
fi

# Check dependencies
check_command "STAR"

# Create output directory
mkdir -p "${output_dir}"

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Starting STAR alignment for: ${sample_id}"
log "INFO" "=========================================="

# Log tool version for reproducibility
star_version="$(STAR --version 2>&1 | head -1)"
log "INFO" "Using STAR ${star_version}"

# Define output prefix (STAR will append suffixes)
output_prefix="${output_dir}/${sample_id}_"

# Build input file string for STAR
# STAR expects space-separated files for R1 and R2
if [[ -n "${input_r2}" ]]; then
    read_files="${input_r1} ${input_r2}"
    log "INFO" "Running in paired-end mode"
else
    read_files="${input_r1}"
    log "INFO" "Running in single-end mode"
fi

# Determine read decompression command based on file extension
read_cmd="cat"
if [[ "${input_r1}" == *.gz ]]; then
    read_cmd="zcat"
fi

log "INFO" "Executing STAR alignment..."
log "INFO" "Genome index: ${genome_index}"
log "INFO" "Threads: ${threads}"

# ==============================================================================
# STAR ALIGNMENT
# ==============================================================================
# The following parameters are optimized for RNA-seq analysis and follow
# Plasmidsaurus and ENCODE recommendations

STAR \
    --runThreadN "${threads}" \
    --genomeDir "${genome_index}" \
    --readFilesIn ${read_files} \
    --readFilesCommand "${read_cmd}" \
    --outFileNamePrefix "${output_prefix}" \
    \
    `# Output format: coordinate-sorted BAM for downstream tools` \
    --outSAMtype BAM SortedByCoordinate \
    \
    `# Filter by splice junctions to remove false positive alignments` \
    `# This is crucial for removing alignments that span non-canonical junctions` \
    --outFilterType BySJout \
    \
    `# Multi-mapping settings: allow up to 20 locations for multi-mappers` \
    `# Multi-mappers are handled downstream by featureCounts with fractional assignment` \
    --outFilterMultimapNmax "${DEFAULT_MULTIMAP_MAX}" \
    \
    `# Mismatch filtering: allow mismatches up to 4% of read length` \
    `# This is permissive enough for SNPs but filters poor alignments` \
    --outFilterMismatchNmax "${DEFAULT_MISMATCH_MAX}" \
    --outFilterMismatchNoverReadLmax "${DEFAULT_MISMATCH_RATIO}" \
    \
    `# Splice junction overhang requirements` \
    `# Unannotated junctions need more support (8bp) to be confident` \
    `# Annotated junctions need minimal support (1bp) as they are known` \
    --alignSJoverhangMin "${DEFAULT_SJ_OVERHANG_UNANNOTATED}" \
    --alignSJDBoverhangMin "${DEFAULT_SJ_OVERHANG_ANNOTATED}" \
    \
    `# Intron size limits to filter implausible alignments` \
    --alignIntronMin "${DEFAULT_INTRON_MIN}" \
    --alignIntronMax "${DEFAULT_INTRON_MAX}" \
    --alignMatesGapMax "${DEFAULT_INTRON_MAX}" \
    \
    `# Output unmapped reads for troubleshooting and QC` \
    `# Useful for detecting contamination or novel sequences` \
    --outReadsUnmapped Fastx \
    \
    `# SAM attributes to include in output` \
    `# NH: number of hits, HI: hit index, AS: alignment score, nM: mismatches` \
    --outSAMattributes NH HI AS nM MD \
    \
    `# Limit BAM sorting memory to prevent memory issues` \
    --limitBAMsortRAM 30000000000

# Check if STAR completed successfully
if [[ $? -ne 0 ]]; then
    die "STAR alignment failed for sample: ${sample_id}"
fi

log "INFO" "STAR alignment completed successfully"

# ==============================================================================
# OUTPUT VERIFICATION
# ==============================================================================
log "INFO" "Verifying outputs..."

bam_file="${output_prefix}Aligned.sortedByCoord.out.bam"
log_file="${output_prefix}Log.final.out"

[[ -s "${bam_file}" ]] || die "Output BAM file is empty or missing: ${bam_file}"
[[ -s "${log_file}" ]] || die "Log.final.out file is empty or missing: ${log_file}"

# ==============================================================================
# ALIGNMENT SUMMARY
# ==============================================================================
log "INFO" "Alignment summary:"

# Parse key metrics from STAR's Log.final.out
if [[ -f "${log_file}" ]]; then
    total_reads=$(grep "Number of input reads" "${log_file}" | cut -f2)
    unique_rate=$(grep "Uniquely mapped reads %" "${log_file}" | cut -f2)
    multi_rate=$(grep "% of reads mapped to multiple loci" "${log_file}" | cut -f2)
    unmapped_rate=$(grep "% of reads unmapped: too short" "${log_file}" | cut -f2)

    log "INFO" "  Input reads:       ${total_reads:-N/A}"
    log "INFO" "  Uniquely mapped:   ${unique_rate:-N/A}"
    log "INFO" "  Multi-mapped:      ${multi_rate:-N/A}"
    log "INFO" "  Unmapped (short):  ${unmapped_rate:-N/A}"
fi

# ==============================================================================
# COMPLETION
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "STAR alignment completed for: ${sample_id}"
log "INFO" "Output BAM: ${bam_file}"
log "INFO" "=========================================="

exit 0
