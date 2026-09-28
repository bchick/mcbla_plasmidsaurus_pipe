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
| 11 | GSEApy | 1.3.1 | Functional enrichment (MSigDB Hallmark) |

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
  - gseapy >= 1.3.1
```

### Reference Files

- Reference genome (FASTA)
- Gene annotation (GTF)
- STAR genome index
- BED12 gene model (RSeQC), and optionally a housekeeping-gene BED12

On the lab server all of these are in `/data/resource` and listed in
`/data/resource/manifest.yaml`: mouse is the Plasmidsaurus portal's Ensembl
114 GRCm39 + ERCC92 build, and human is GENCODE v44.

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

**On the Salk lab server** the references are already in `/data/resource`.
`pixi run init` asks which genome you are using and writes a project config
that points at them:

```bash
pixi install
pixi run init --list                                   # hg38 (GENCODE v44) or mm39 (Ensembl 114 + ERCC92)
pixi run init --dir /data/<user>/<project> --genome mm39
# edit <project>/samples.tsv and edger.contrasts in <project>/config.yaml, then
pixi run bash run_pipeline.sh -i /path/to/fastq -o /data/<user>/<project>/results \
    -c /data/<user>/<project>/config.yaml -m /data/<user>/<project>/samples.tsv --dry-run
```

Agents running the pipeline for someone should follow [AGENTS.md](AGENTS.md).

**Elsewhere:** copy `config/config.template.yaml`, fill in your reference
paths (FASTA, GTF, STAR index, BED12 gene model), and pass it with `-c`.

## Usage

```bash
pixi run bash run_pipeline.sh -i <input_dir> -o <output_dir> [-g genome | -c config] [OPTIONS]

  -i, --input       directory with the FASTQ (or BAM) files
  -o, --output      output directory
  -g, --genome      hg38 (config/config.human.yaml) or mm39 (config/config.mouse.yaml);
                    mm10 is a deprecated alias of mm39
  -c, --config      YAML config (instead of -g), e.g. the one `pixi run init` wrote
  -m, --metadata    sample sheet (sample_id, condition, replicate); needed for steps 8-10
  -y, --type        fastq (default) or bam (starts at step 5)
  -s, --start-step / -e, --end-step   steps 1-10
  -t, --threads     threads (default 8)
  -d, --dry-run     print the commands only
```

See `./run_pipeline.sh -h` for the full help.

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
