# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/claude-code) when working with this repository.

## Project Overview

This is a bioinformatics pipeline for processing and analyzing RNA-seq data from Plasmidsaurus sequencing services. This is **not** for plasmid sequencing/verification - it processes RNA-seq data.

## Tech Stack

- **Primary language**: Bash/Shell scripts
- **Domain**: Bioinformatics / RNA-seq analysis

## Pipeline Steps

1. **FastQ generation & demux**: BCL Convert v4.3.6, fqtk v0.3.1
2. **Read filtering**: FastP v0.24.0
   - Poly-X tail trimming
   - 3' quality-based tail trimming
   - Min Phred quality: 15
   - Min read length: 50 bp
3. **Alignment**: STAR v2.7.11
   - Non-canonical splice junction removal
   - Output unmapped reads
4. **BAM sorting**: samtools v1.22.1 (coordinate sort)
5. **UMI deduplication**: UMICollapse v1.1.0 (PCR + optical duplicate removal)
6. **Mapping QC**: RSeQC v5.0.4, Qualimap v2.3
   - Alignment quality metrics
   - Strand specificity
   - Read distribution across genomic features
7. **QC reporting**: MultiQC v1.32
8. **Gene quantification**: featureCounts (subread v2.1.1)
   - Strand-specific counting
   - Multi-mapping read fractional assignment
   - Features: exons and 3' UTR
   - Grouped by gene_id
   - Annotated with gene biotype and GTF metadata
9. **Sample correlations**: TMM normalization, Pearson correlation (for heatmap/PCA)
10. **Differential expression**: edgeR v4.0.16
    - Low-expression filtering: `edgeR::filterByExpr` (default values)
11. **Functional enrichment**: GSEApy v0.12
    - MSigDB Hallmark gene set
    - Human and mouse samples

## Common Commands

```bash
# Make scripts executable
chmod +x *.sh

# Run the main pipeline (update as developed)
# ./run_pipeline.sh <input_dir> <output_dir>
```

## Project Structure

```
├── scripts/          # Individual analysis scripts
├── data/             # Input data (not tracked in git)
├── results/          # Output results (not tracked in git)
├── config/           # Configuration files
└── logs/             # Pipeline logs
```

## Code Quality Standards

**All code must be written to professional, production-ready standards.** The goal is that even experienced bioinformaticians would not need to redo any steps.

### Documentation Requirements

Every script MUST include:

1. **File header block** with:
   - Script name and brief description
   - Author and date
   - Usage examples
   - Input/output specifications
   - Dependencies and version requirements

2. **Section comments** that divide the script into logical blocks:
   ```bash
   # ==============================================================================
   # SECTION NAME
   # ==============================================================================
   # Brief description of what this section accomplishes
   ```

3. **Inline comments** explaining:
   - Why a particular approach was chosen (not just what it does)
   - Non-obvious parameter choices and their biological rationale
   - Edge cases being handled
   - Any assumptions being made

4. **Function documentation**:
   ```bash
   #######################################
   # Brief description of function purpose.
   #
   # Globals:
   #   VARIABLE_NAME - description of global used
   # Arguments:
   #   $1 - description of first argument
   #   $2 - description of second argument
   # Outputs:
   #   Writes to stdout/file description
   # Returns:
   #   0 on success, non-zero on error
   #######################################
   function_name() {
       ...
   }
   ```

### Bash Script Standards

```bash
#!/usr/bin/env bash
#
# script_name.sh - Brief description of script purpose
#
# Description:
#   Detailed description of what this script does, why it exists,
#   and how it fits into the overall pipeline.
#
# Usage:
#   ./script_name.sh -i <input_file> -o <output_dir> [-t threads]
#
# Arguments:
#   -i, --input     Path to input FASTQ file (required)
#   -o, --output    Path to output directory (required)
#   -t, --threads   Number of threads to use (default: 4)
#   -h, --help      Display this help message
#
# Dependencies:
#   - tool_name >= version (purpose)
#   - another_tool >= version (purpose)
#
# Example:
#   ./script_name.sh -i sample.fastq.gz -o results/ -t 8
#
# Author: Your Name
# Date: YYYY-MM-DD
# Version: 1.0.0

# Exit immediately on error, undefined variable, or pipe failure
set -euo pipefail

# Enable debug mode if DEBUG environment variable is set
[[ "${DEBUG:-}" == "true" ]] && set -x
```

### Variable and Function Naming

- **Local variables**: `lowercase_with_underscores`
- **Environment/global variables**: `UPPERCASE_WITH_UNDERSCORES`
- **Constants**: `readonly CONSTANT_NAME="value"`
- **Functions**: `lowercase_with_underscores()`
- **Temporary files**: Use `mktemp` with descriptive prefixes

### Error Handling

```bash
# Always provide informative error messages
die() {
    echo "ERROR: $*" >&2
    exit 1
}

# Validate inputs early
[[ -f "${input_file}" ]] || die "Input file not found: ${input_file}"
[[ -d "${output_dir}" ]] || mkdir -p "${output_dir}"

# Use trap for cleanup
trap 'cleanup' EXIT ERR
```

### Logging Standards

```bash
# Timestamp function for consistent log formatting
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

log "INFO: Starting alignment step..."
log "WARN: Low read count detected in sample ${sample_id}"
log "ERROR: Alignment failed for ${input_file}"
```

### Command Documentation

When calling bioinformatics tools, document key parameters:

```bash
# Run STAR alignment
# --outFilterType BySJout: Filter by splice junctions for cleaner output
# --outFilterMultimapNmax 20: Allow up to 20 multi-mapping locations
# --alignSJoverhangMin 8: Minimum overhang for unannotated junctions
# --outSAMtype BAM SortedByCoordinate: Output coordinate-sorted BAM
STAR \
    --runThreadN "${threads}" \
    --genomeDir "${genome_index}" \
    --readFilesIn "${input_r1}" "${input_r2}" \
    --readFilesCommand zcat \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${output_prefix}"
```

## Input/Output

- **Input**: FASTQ files from Plasmidsaurus RNA-seq
- **Output**: Alignment files, count matrices, QC reports, differential expression results, enrichment analysis

## Important Notes

- Always validate input file existence before processing
- Log all command versions for reproducibility
- Handle edge cases (empty files, failed alignments)
- Write code as if the next person to read it has no context - because they won't
