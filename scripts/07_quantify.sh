#!/usr/bin/env bash
#
# 07_quantify.sh - Quantify gene expression using featureCounts
#
# Description:
#   This script performs gene-level expression quantification using featureCounts
#   from the Subread package. By default it counts reads over whole gene bodies
#   (GTF "gene" features), grouped by gene ID, which is what the Plasmidsaurus
#   portal does. Multi-mapping reads are counted fractionally (1/NH each),
#   also as the portal does.
#
# Usage:
#   ./07_quantify.sh -i <input_bams> -g <gtf_file> -o <output_dir> [OPTIONS]
#
# Arguments:
#   -i, --input         Input BAM file(s), comma-separated or directory (required)
#   -g, --gtf           Path to GTF annotation file (required)
#   -o, --output-dir    Output directory for count matrices (required)
#   -s, --strand        Strandedness: 0=unstranded, 1=forward, 2=reverse (default: 1)
#   -t, --threads       Number of threads (default: 8)
#   -f, --feature       Feature type to count (default: gene)
#   -q, --min-mapq      Minimum mapping quality (default: 0)
#   -a, --attribute     Attribute for grouping (default: gene_id)
#   -h, --help          Display this help message
#
# Dependencies:
#   - featureCounts (Subread >= 2.1.1)
#
# Output:
#   - gene_counts.txt           Raw count matrix
#   - gene_counts.txt.summary   Counting summary statistics
#   - gene_counts_annotated.txt Count matrix with gene annotations
#
# Example:
#   ./07_quantify.sh -i results/04_dedup/ -g annotation.gtf -o results/07_counts/ -s 1
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
# Plasmidsaurus reads are sense-stranded (forward, -s 1); verified on run
# LJQQSK, where -s 1 reproduces the portal count matrix exactly.
readonly DEFAULT_STRAND=1
# Count over whole gene bodies (introns included), as the portal does. With
# "exon", intronic reads from pre-mRNA are dropped and counts fall short of
# the portal's.
readonly DEFAULT_FEATURE="gene"
readonly DEFAULT_ATTRIBUTE="gene_id"
# STAR assigns MAPQ 255 to unique reads and 0/1/3 to multi-mappers, so any
# threshold above 0 silently discards every multi-mapper and makes -M --fraction
# a no-op. Plasmidsaurus counts multi-mappers fractionally, so keep them all.
readonly DEFAULT_MIN_MAPQ=0        # Minimum mapping quality

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
Usage: ${SCRIPT_NAME} -i <input_bams> -g <gtf_file> -o <output_dir> [OPTIONS]

Quantify gene expression using featureCounts.

Required Arguments:
  -i, --input         Input BAM file(s): comma-separated list, or directory path
                      (will find all .bam files in directory)
  -g, --gtf           Path to GTF annotation file
  -o, --output-dir    Output directory for count matrices

Optional Arguments:
  -s, --strand        Strandedness: 0=unstranded, 1=forward, 2=reverse (default: ${DEFAULT_STRAND})
  -t, --threads       Number of threads (default: ${DEFAULT_THREADS})
  -f, --feature       Feature type to count (default: ${DEFAULT_FEATURE})
  -a, --attribute     Attribute for grouping features (default: ${DEFAULT_ATTRIBUTE})
  -q, --min-mapq      Minimum mapping quality (default: ${DEFAULT_MIN_MAPQ})
  -h, --help          Display this help message

Strandedness Options:
  0 - Unstranded (count reads regardless of strand)
  1 - Forward/sense stranded (read strand matches gene strand)
      Plasmidsaurus RNA-seq libraries (default)
  2 - Reverse/antisense stranded (read strand opposite to gene strand)
      Most common for Illumina dUTP-based stranded library prep

Examples:
  # Count reads from all BAMs in a directory (forward stranded)
  ${SCRIPT_NAME} -i results/04_dedup/ -g annotation.gtf -o results/07_counts/ -s 1

  # Count from specific BAM files
  ${SCRIPT_NAME} -i "sample1.bam,sample2.bam" -g annotation.gtf -o results/07_counts/

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

input_bams=""
gtf_file=""
output_dir=""
strand="${DEFAULT_STRAND}"
threads="${DEFAULT_THREADS}"
feature="${DEFAULT_FEATURE}"
attribute="${DEFAULT_ATTRIBUTE}"
min_mapq="${DEFAULT_MIN_MAPQ}"
paired_end="auto"  # auto, yes, no - auto-detect by default

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input_bams="$2"
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
        -s|--strand)
            strand="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -f|--feature)
            feature="$2"
            shift 2
            ;;
        -a|--attribute)
            attribute="$2"
            shift 2
            ;;
        -q|--min-mapq)
            min_mapq="$2"
            shift 2
            ;;
        -p|--paired)
            paired_end="$2"
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

[[ -z "${input_bams}" ]] && die "Input BAM file(s) required (-i)"
[[ -z "${gtf_file}" ]] && die "GTF annotation file required (-g)"
[[ -z "${output_dir}" ]] && die "Output directory required (-o)"
[[ -f "${gtf_file}" ]] || die "GTF file not found: ${gtf_file}"

# Validate strand parameter
case "${strand}" in
    0|1|2)
        ;;
    *)
        die "Invalid strand value: ${strand}. Must be 0, 1, or 2."
        ;;
esac

check_command "featureCounts"
mkdir -p "${output_dir}"

# ==============================================================================
# RESOLVE INPUT BAM FILES
# ==============================================================================
log "INFO" "Resolving input BAM files..."

bam_files=()

# Check if input is a directory
if [[ -d "${input_bams}" ]]; then
    # Find all BAM files in directory
    while IFS= read -r -d '' bam; do
        bam_files+=("${bam}")
    done < <(find "${input_bams}" -name "*.bam" -type f -print0 | sort -z)

    [[ ${#bam_files[@]} -eq 0 ]] && die "No BAM files found in directory: ${input_bams}"
else
    # Parse comma-separated list
    IFS=',' read -ra bam_array <<< "${input_bams}"
    for bam in "${bam_array[@]}"; do
        bam=$(echo "${bam}" | xargs)  # Trim whitespace
        [[ -f "${bam}" ]] || die "BAM file not found: ${bam}"
        bam_files+=("${bam}")
    done
fi

log "INFO" "Found ${#bam_files[@]} BAM file(s)"

# ==============================================================================
# DETECT PAIRED-END STATUS
# ==============================================================================
# Auto-detect paired-end status from the first BAM file if not specified
if [[ "${paired_end}" == "auto" ]]; then
    first_bam="${bam_files[0]}"
    log "INFO" "Auto-detecting paired-end status from: $(basename "${first_bam}")"

    # Use samtools to check if reads are paired
    # Check the FLAG field of the first 1000 reads for paired-end flag (0x1)
    # Use subshell to avoid SIGPIPE issues with head
    paired_count=$(set +o pipefail; samtools view "${first_bam}" 2>/dev/null | head -1000 | awk '$2 % 2 == 1 {count++} END {print count+0}')

    if [[ "${paired_count}" -gt 500 ]]; then
        paired_end="yes"
        log "INFO" "Detected: Paired-end reads (${paired_count}/1000 reads are paired)"
    else
        paired_end="no"
        log "INFO" "Detected: Single-end reads (${paired_count}/1000 reads are paired)"
    fi
fi

# ==============================================================================
# MAIN PROCESSING
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Starting gene quantification"
log "INFO" "=========================================="

# Log tool version
# featureCounts -v prints a leading blank line, so match the version line itself
fc_version="$(featureCounts -v 2>&1 | grep -m1 'v[0-9]')"
log "INFO" "Using ${fc_version}"

# Log parameters
log "INFO" "Parameters:"
log "INFO" "  GTF file:     ${gtf_file}"
log "INFO" "  Feature type: ${feature}"
log "INFO" "  Attribute:    ${attribute}"
log "INFO" "  Strandedness: ${strand} (0=unstranded, 1=forward, 2=reverse)"
log "INFO" "  Paired-end:   ${paired_end}"
log "INFO" "  Min MAPQ:     ${min_mapq}"
log "INFO" "  Threads:      ${threads}"

output_counts="${output_dir}/gene_counts.txt"

# ==============================================================================
# RUN FEATURECOUNTS
# ==============================================================================
log "INFO" "Running featureCounts..."

# featureCounts parameters explained:
# -a GTF: Annotation file
# -o output: Output file path
# -t feature: Feature type (gene = whole gene body, the portal method)
# -g attribute: Attribute for grouping (gene_id)
# -s strand: Strandedness setting
# -T threads: Number of threads
# -p: Count fragments instead of reads (for paired-end ONLY)
# -B: Only count read pairs with both ends mapped (paired-end ONLY)
# -C: Don't count pairs mapping to different chromosomes (paired-end ONLY)
# -Q min_mapq: Minimum mapping quality
# -M: Count multi-mapping reads
# --fraction: Assign fractional counts for multi-mappers
#             (e.g., if read maps to 4 genes, each gets 0.25)
# (no -O): Reads overlapping more than one gene are left unassigned as
#          ambiguous (featureCounts default) rather than split between genes
# --extraAttributes: Include additional attributes from GTF

# Ensembl GTFs call the biotype attribute "gene_biotype"; GENCODE calls it
# "gene_type". Detect which one this GTF uses so the column is populated.
# Use a subshell with pipefail off: grep -m1 closes the pipe early (SIGPIPE).
if (set +o pipefail; zcat -f "${gtf_file}" | grep -m1 -v '^#' | grep -q 'gene_type "'); then
    biotype_attr="gene_type"
else
    biotype_attr="gene_biotype"
fi
log "INFO" "  Biotype attr: ${biotype_attr}"

# Build featureCounts command with conditional paired-end flags
fc_cmd=(featureCounts
    -a "${gtf_file}"
    -o "${output_counts}"
    -t "${feature}"
    -g "${attribute}"
    -s "${strand}"
    -T "${threads}"
    -Q "${min_mapq}"
    -M
    --fraction
    --extraAttributes "gene_name,${biotype_attr}"
)

# Add paired-end specific flags only if data is paired-end
if [[ "${paired_end}" == "yes" ]]; then
    fc_cmd+=(-p -B -C)
    # Since subread 2.0.2, -p alone counts each mate separately; fragments are
    # only counted when --countReadPairs is also given. Older versions reject it.
    fc_ver_num="$(sed -E 's/.*v([0-9.]+).*/\1/' <<< "${fc_version}")"
    if [[ "$(printf '%s\n' 2.0.2 "${fc_ver_num}" | sort -V | head -1)" == "2.0.2" ]]; then
        fc_cmd+=(--countReadPairs)
    fi
    log "INFO" "Using paired-end mode (counting fragments)"
else
    log "INFO" "Using single-end mode"
fi

# Add BAM files and run
"${fc_cmd[@]}" "${bam_files[@]}"

log "INFO" "featureCounts completed"

# ==============================================================================
# POST-PROCESSING
# ==============================================================================
log "INFO" "Post-processing count matrix..."

# Downstream scripts expect the biotype column to be named gene_biotype
if [[ "${biotype_attr}" != "gene_biotype" ]]; then
    sed -i "2s/\t${biotype_attr}\t/\tgene_biotype\t/" "${output_counts}"
fi

# Create a cleaner version of the counts file
# Remove the first comment line and simplify column headers
annotated_counts="${output_dir}/gene_counts_annotated.txt"

# Extract header and clean up BAM file paths to sample names
head -2 "${output_counts}" | tail -1 | \
    sed 's|[^\t]*/||g' | \
    sed 's|\.dedup\.bam||g' | \
    sed 's|\.sorted\.bam||g' | \
    sed 's|\.bam||g' \
    > "${annotated_counts}"

# Add data rows (skip header lines)
tail -n +3 "${output_counts}" >> "${annotated_counts}"

log "INFO" "Annotated count matrix: ${annotated_counts}"

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
log "INFO" "Generating summary statistics..."

summary_file="${output_counts}.summary"
if [[ -f "${summary_file}" ]]; then
    log "INFO" "Assignment summary:"

    # Parse and display key statistics
    while IFS=$'\t' read -r status count rest; do
        if [[ "${status}" == "Assigned" ]]; then
            log "INFO" "  Assigned:              ${count}"
        elif [[ "${status}" == "Unassigned_Unmapped" ]]; then
            log "INFO" "  Unassigned (unmapped): ${count}"
        elif [[ "${status}" == "Unassigned_NoFeatures" ]]; then
            log "INFO" "  Unassigned (no feat):  ${count}"
        elif [[ "${status}" == "Unassigned_Ambiguity" ]]; then
            log "INFO" "  Unassigned (ambig):    ${count}"
        fi
    done < "${summary_file}"
fi

# Count genes with expression
genes_detected=$(tail -n +3 "${output_counts}" | awk -F'\t' '
    {
        # Sum counts across all samples (columns 7+)
        sum = 0
        for (i = 7; i <= NF; i++) sum += $i
        if (sum > 0) detected++
    }
    END { print detected }
')
total_genes=$(tail -n +3 "${output_counts}" | wc -l)

log "INFO" "Gene detection:"
log "INFO" "  Total genes in GTF:    ${total_genes}"
log "INFO" "  Genes with reads:      ${genes_detected}"

# ==============================================================================
# COMPLETION
# ==============================================================================
log "INFO" "=========================================="
log "INFO" "Gene quantification completed"
log "INFO" "Count matrix:    ${output_counts}"
log "INFO" "Summary:         ${summary_file}"
log "INFO" "=========================================="

exit 0
