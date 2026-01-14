# Plasmidsaurus RNA-seq Analysis Pipeline

A robust, production-ready bioinformatics pipeline for processing and analyzing RNA-seq data from Plasmidsaurus sequencing services.

> **Note**: This pipeline is for RNA-seq analysis, not plasmid sequencing/verification.

## Overview

This pipeline implements a comprehensive RNA-seq workflow from raw FASTQ files through differential expression and functional enrichment analysis. The workflow is designed to be reproducible, well-documented, and suitable for production use.

## Pipeline Steps

| Step | Tool | Version | Description |
|------|------|---------|-------------|
| 1 | BCL Convert / fqtk | 4.3.6 / 0.3.1 | FastQ generation and demultiplexing |
| 2 | FastP | 0.24.0 | Read filtering and QC |
| 3 | STAR | 2.7.11 | Alignment to reference genome |
| 4 | samtools | 1.22.1 | BAM coordinate sorting |
| 5 | UMICollapse | 1.1.0 | UMI-based deduplication |
| 6 | RSeQC / Qualimap | 5.0.4 / 2.3 | Mapping quality control |
| 7 | MultiQC | 1.32 | Comprehensive QC reporting |
| 8 | featureCounts | 2.1.1 | Gene expression quantification |
| 9 | edgeR | 4.0.16 | TMM normalization and sample correlations |
| 10 | edgeR | 4.0.16 | Differential expression analysis |
| 11 | GSEApy | 0.12 | Functional enrichment (MSigDB Hallmark) |

## Requirements

### Software Dependencies

```
bcl-convert >= 4.3.6
fqtk >= 0.3.1
fastp >= 0.24.0
STAR >= 2.7.11
samtools >= 1.22.1
UMICollapse >= 1.1.0
rseqc >= 5.0.4
qualimap >= 2.3
multiqc >= 1.32
subread >= 2.1.1 (featureCounts)
R >= 4.0
  - edgeR >= 4.0.16
  - DESeq2 (optional)
Python >= 3.8
  - gseapy >= 0.12
```

### Reference Files

- Reference genome (FASTA)
- Gene annotation (GTF)
- STAR genome index

## Project Structure

```
├── scripts/          # Pipeline scripts
├── config/           # Configuration files
├── data/             # Input data (not tracked)
├── results/          # Output results (not tracked)
├── logs/             # Execution logs (not tracked)
├── CLAUDE.md         # Claude Code guidance
└── README.md         # This file
```

## Quick Start

```bash
# 1. Clone the repository
git clone <repository-url>
cd mcbla_plasmidsaurus_pipe

# 2. Configure the pipeline
cp config/config.template.yaml config/config.yaml
# Edit config/config.yaml with your paths and parameters

# 3. Run the pipeline
./run_pipeline.sh -i /path/to/fastq -o results/ -c config/config.yaml
```

## Usage

```bash
./run_pipeline.sh [OPTIONS]

Options:
  -i, --input DIR       Input directory containing FASTQ files (required)
  -o, --output DIR      Output directory for results (required)
  -c, --config FILE     Configuration file (default: config/config.yaml)
  -t, --threads INT     Number of threads (default: 8)
  -s, --step STEP       Start from specific step (default: 1)
  -h, --help            Display help message

Examples:
  # Run full pipeline
  ./run_pipeline.sh -i data/fastq/ -o results/ -t 16

  # Resume from alignment step
  ./run_pipeline.sh -i data/fastq/ -o results/ -s 3
```

## Output Structure

```
results/
├── 01_filtered/           # FastP filtered reads
├── 02_aligned/            # STAR alignments
├── 03_sorted/             # Coordinate-sorted BAMs
├── 04_dedup/              # Deduplicated BAMs
├── 05_qc/                 # RSeQC and Qualimap output
├── 06_multiqc/            # MultiQC report
├── 07_counts/             # featureCounts matrices
├── 08_normalization/      # TMM-normalized counts
├── 09_de/                 # Differential expression results
└── 10_enrichment/         # GSEA results
```

## Configuration

See `config/config.template.yaml` for all available options including:

- Reference genome and annotation paths
- Tool-specific parameters
- Sample metadata
- Contrast definitions for DE analysis

## License

[Add license information]

## Authors

[Add author information]

## Acknowledgments

Pipeline design based on Plasmidsaurus bioinformatics workflows.
