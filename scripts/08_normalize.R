#!/usr/bin/env Rscript
#
# 08_normalize.R - Normalize counts and generate sample correlation metrics
#
# Description:
#   This script performs TMM (Trimmed Mean of M-values) normalization on raw
#   gene counts using edgeR. TMM normalization accounts for composition bias
#   between samples, which can occur when a small number of highly expressed
#   genes consume a large proportion of sequencing reads. The script also
#   generates sample-sample correlations for quality assessment via heatmaps
#   and PCA plots.
#
# Usage:
#   Rscript 08_normalize.R --counts <count_matrix> --metadata <sample_info> \
#       --output <output_dir>
#
# Arguments:
#   --counts      Path to gene count matrix from featureCounts (required)
#   --metadata    Path to sample metadata TSV file (required)
#   --output      Output directory for normalized data (required)
#   --min-count   Minimum count threshold for filtering (default: 10)
#   --help        Display this help message
#
# Dependencies:
#   - R >= 4.0
#   - edgeR >= 4.0.16
#   - ggplot2
#   - pheatmap
#   - RColorBrewer
#
# Output:
#   - normalized_counts_tmm.txt    TMM-normalized counts (CPM)
#   - normalized_counts_log2.txt   Log2-transformed normalized counts
#   - sample_correlation.txt       Pearson correlation matrix
#   - plots/pca_plot.pdf           PCA plot
#   - plots/correlation_heatmap.pdf Sample correlation heatmap
#   - plots/sample_distances.pdf   Sample distance heatmap
#
# Example:
#   Rscript 08_normalize.R \
#       --counts results/07_counts/gene_counts.txt \
#       --metadata config/samples.tsv \
#       --output results/08_normalization/
#
# Author: Plasmidsaurus RNA-seq Pipeline
# Date: 2024
# Version: 1.0.0

# ==============================================================================
# SETUP AND LIBRARIES
# ==============================================================================

# Suppress package startup messages for cleaner output
suppressPackageStartupMessages({
    library(optparse)
    library(edgeR)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

# Set seed for reproducibility (affects some visualization elements)
set.seed(42)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Log message with timestamp
#'
#' @param level Log level (INFO, WARN, ERROR)
#' @param message Message to log
log_message <- function(level, message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", timestamp, level, message), file = stderr())
}

#' Stop execution with error message
#'
#' @param message Error message
die <- function(message) {
    log_message("ERROR", message)
    quit(status = 1)
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

option_list <- list(
    make_option(c("--counts"), type = "character", default = NULL,
                help = "Path to gene count matrix from featureCounts [required]"),
    make_option(c("--metadata"), type = "character", default = NULL,
                help = "Path to sample metadata TSV file [required]"),
    make_option(c("--output"), type = "character", default = NULL,
                help = "Output directory for normalized data [required]"),
    make_option(c("--min-count"), type = "integer", default = 10,
                help = "Minimum count threshold for filtering [default: %default]")
)

# Parse arguments
opt_parser <- OptionParser(
    option_list = option_list,
    description = "Normalize gene counts using TMM and generate sample correlations"
)
opt <- parse_args(opt_parser)

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_message("INFO", "==========================================")
log_message("INFO", "Starting count normalization")
log_message("INFO", "==========================================")

# Check required arguments
if (is.null(opt$counts)) die("Count matrix file is required (--counts)")
if (is.null(opt$metadata)) die("Sample metadata file is required (--metadata)")
if (is.null(opt$output)) die("Output directory is required (--output)")

# Check input files exist
if (!file.exists(opt$counts)) die(paste("Count matrix not found:", opt$counts))
if (!file.exists(opt$metadata)) die(paste("Metadata file not found:", opt$metadata))

# Create output directories
dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(opt$output, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

log_message("INFO", paste("Count matrix:", opt$counts))
log_message("INFO", paste("Metadata:", opt$metadata))
log_message("INFO", paste("Output directory:", opt$output))

# ==============================================================================
# LOAD DATA
# ==============================================================================

log_message("INFO", "Loading count matrix...")

# Read featureCounts output
# featureCounts format: Geneid, Chr, Start, End, Strand, Length, Sample1, Sample2, ...
# We skip the first row (comment line starting with #) and use row 2 as header
counts_raw <- read.table(
    opt$counts,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    comment.char = "#",
    check.names = FALSE
)

# Extract count columns (remove annotation columns: Chr, Start, End, Strand, Length)
# Also handle extraAttributes columns (gene_name, gene_biotype) if present
annotation_cols <- c("Chr", "Start", "End", "Strand", "Length", "gene_name", "gene_biotype")
annotation_cols <- annotation_cols[annotation_cols %in% colnames(counts_raw)]

if (length(annotation_cols) > 0) {
    gene_annotations <- counts_raw[, annotation_cols, drop = FALSE]
    counts <- counts_raw[, !colnames(counts_raw) %in% annotation_cols, drop = FALSE]
} else {
    gene_annotations <- NULL
    counts <- counts_raw
}

# Clean sample names (remove path and extension artifacts)
colnames(counts) <- gsub(".*\\/", "", colnames(counts))
colnames(counts) <- gsub("_dedup-mapped-reads\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("\\.dedup\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("\\.sorted\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("\\.bam$", "", colnames(counts))

log_message("INFO", paste("Loaded counts for", nrow(counts), "genes and", ncol(counts), "samples"))

# Load sample metadata
log_message("INFO", "Loading sample metadata...")
metadata <- read.table(
    opt$metadata,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Verify sample IDs match between counts and metadata
if (!"sample_id" %in% colnames(metadata)) {
    die("Metadata must contain 'sample_id' column")
}

# Match metadata to count matrix column order
metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
if (any(is.na(metadata$sample_id))) {
    missing <- colnames(counts)[!colnames(counts) %in% metadata$sample_id]
    die(paste("Samples in count matrix not found in metadata:", paste(missing, collapse = ", ")))
}

log_message("INFO", paste("Metadata loaded for", nrow(metadata), "samples"))

# ==============================================================================
# CREATE DGELIST OBJECT
# ==============================================================================

log_message("INFO", "Creating DGEList object...")

# DGEList is edgeR's main data structure
# It stores counts, sample information, and normalization factors
dge <- DGEList(
    counts = as.matrix(counts),
    samples = metadata,
    genes = gene_annotations
)

log_message("INFO", paste("Initial genes:", nrow(dge)))
log_message("INFO", paste("Initial samples:", ncol(dge)))

# ==============================================================================
# FILTER LOW-EXPRESSED GENES
# ==============================================================================

log_message("INFO", "Filtering low-expressed genes...")

# filterByExpr uses edgeR's recommended filtering strategy
# It keeps genes with sufficient counts in enough samples
# The filtering is based on:
# 1. Minimum count threshold (min.count)
# 2. Minimum total count across all samples (min.total.count)
# 3. CPM threshold derived from smallest library size
# 4. Proportion of samples that must meet the threshold

# Determine group variable for filtering (if condition column exists)
group <- NULL
if ("condition" %in% colnames(metadata)) {
    group <- factor(metadata$condition)
    log_message("INFO", paste("Using 'condition' for group-aware filtering:",
                              paste(levels(group), collapse = ", ")))
}

# Apply filtering
keep <- filterByExpr(
    dge,
    group = group,
    min.count = opt$`min-count`
)

dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

genes_removed <- nrow(dge) - nrow(dge_filtered)
log_message("INFO", paste("Genes removed (low expression):", genes_removed))
log_message("INFO", paste("Genes retained:", nrow(dge_filtered)))

# ==============================================================================
# TMM NORMALIZATION
# ==============================================================================

log_message("INFO", "Performing TMM normalization...")

# calcNormFactors calculates normalization factors using the TMM method
# TMM (Trimmed Mean of M-values) method:
# 1. Selects a reference sample (usually the one with upper quartile closest to mean)
# 2. Calculates M-values (log fold changes) between each sample and reference
# 3. Trims extreme M-values (default: 30% from each tail)
# 4. Calculates weighted mean of remaining M-values as the normalization factor
#
# This approach is robust to:
# - Highly expressed genes that dominate library composition
# - Genes with very low or zero counts
# - Asymmetric differential expression (more up than down, or vice versa)

dge_normalized <- calcNormFactors(dge_filtered, method = "TMM")

# Log normalization factors
log_message("INFO", "TMM normalization factors:")
for (i in seq_len(ncol(dge_normalized))) {
    log_message("INFO", sprintf("  %s: %.3f",
                                colnames(dge_normalized)[i],
                                dge_normalized$samples$norm.factors[i]))
}

# ==============================================================================
# CALCULATE NORMALIZED COUNTS
# ==============================================================================

log_message("INFO", "Calculating normalized expression values...")

# CPM (Counts Per Million) normalizes for library size
# With norm.factors applied, this gives TMM-normalized CPM
cpm_normalized <- cpm(dge_normalized, normalized.lib.sizes = TRUE)

# Log2 transformation for downstream visualization
# Add small offset to avoid log(0)
log2_cpm <- log2(cpm_normalized + 1)

# ==============================================================================
# SAMPLE CORRELATIONS
# ==============================================================================

log_message("INFO", "Calculating sample correlations...")

# Calculate Pearson correlation matrix on log2-transformed normalized counts
# Pearson correlation measures linear relationship between samples
# High correlation (>0.9) between biological replicates indicates good reproducibility
sample_cor <- cor(log2_cpm, method = "pearson")

log_message("INFO", "Sample correlation summary:")
log_message("INFO", sprintf("  Min correlation: %.3f", min(sample_cor[lower.tri(sample_cor)])))
log_message("INFO", sprintf("  Max correlation: %.3f", max(sample_cor[lower.tri(sample_cor)])))
log_message("INFO", sprintf("  Mean correlation: %.3f", mean(sample_cor[lower.tri(sample_cor)])))

# ==============================================================================
# SAVE OUTPUT FILES
# ==============================================================================

log_message("INFO", "Saving output files...")

# Save TMM-normalized CPM
cpm_output <- file.path(opt$output, "normalized_counts_tmm.txt")
cpm_df <- data.frame(gene_id = rownames(cpm_normalized), cpm_normalized, check.names = FALSE)
write.table(cpm_df, cpm_output, sep = "\t", quote = FALSE, row.names = FALSE)
log_message("INFO", paste("TMM-normalized CPM saved:", cpm_output))

# Save log2-transformed counts
log2_output <- file.path(opt$output, "normalized_counts_log2.txt")
log2_df <- data.frame(gene_id = rownames(log2_cpm), log2_cpm, check.names = FALSE)
write.table(log2_df, log2_output, sep = "\t", quote = FALSE, row.names = FALSE)
log_message("INFO", paste("Log2 CPM saved:", log2_output))

# Save correlation matrix
cor_output <- file.path(opt$output, "sample_correlation.txt")
cor_df <- data.frame(sample_id = rownames(sample_cor), sample_cor, check.names = FALSE)
write.table(cor_df, cor_output, sep = "\t", quote = FALSE, row.names = FALSE)
log_message("INFO", paste("Correlation matrix saved:", cor_output))

# Save DGEList object for downstream DE analysis
rds_output <- file.path(opt$output, "dge_normalized.rds")
saveRDS(dge_normalized, rds_output)
log_message("INFO", paste("DGEList object saved:", rds_output))

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

log_message("INFO", "Generating visualization plots...")

# Define color palette for conditions (if available)
if ("condition" %in% colnames(metadata)) {
    n_conditions <- length(unique(metadata$condition))
    condition_colors <- brewer.pal(max(3, n_conditions), "Set1")[1:n_conditions]
    names(condition_colors) <- unique(metadata$condition)
    annotation_col <- data.frame(
        Condition = factor(metadata$condition),
        row.names = metadata$sample_id
    )
    annotation_colors <- list(Condition = condition_colors)
} else {
    annotation_col <- NULL
    annotation_colors <- NULL
}

# ------------------------------------------------------------------
# PCA Plot
# ------------------------------------------------------------------
log_message("INFO", "Creating PCA plot...")

# Perform PCA on log2-transformed normalized counts
# We transpose because prcomp expects samples in rows
pca_result <- prcomp(t(log2_cpm), scale. = TRUE, center = TRUE)

# Calculate variance explained by each PC
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Create PCA data frame for plotting
pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample_id = rownames(pca_result$x)
)

# Add metadata if available
if ("condition" %in% colnames(metadata)) {
    pca_df$condition <- metadata$condition[match(pca_df$sample_id, metadata$sample_id)]
}

# Generate PCA plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = if ("condition" %in% colnames(pca_df)) condition else NULL),
               size = 4, alpha = 0.8) +
    geom_text(aes(label = sample_id), vjust = -1, size = 3) +
    labs(
        title = "PCA of Normalized Gene Expression",
        x = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
        y = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
        color = "Condition"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom"
    )

if ("condition" %in% colnames(pca_df)) {
    pca_plot <- pca_plot + scale_color_manual(values = condition_colors)
}

# Save PCA plot
pca_file <- file.path(plots_dir, "pca_plot.pdf")
ggsave(pca_file, pca_plot, width = 10, height = 8)
log_message("INFO", paste("PCA plot saved:", pca_file))

# ------------------------------------------------------------------
# Correlation Heatmap
# ------------------------------------------------------------------
log_message("INFO", "Creating correlation heatmap...")

# Generate correlation heatmap using pheatmap
cor_heatmap_file <- file.path(plots_dir, "correlation_heatmap.pdf")

pdf(cor_heatmap_file, width = 10, height = 8)
pheatmap(
    sample_cor,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(brewer.pal(9, "Blues"))(100),
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 8,
    main = "Sample-Sample Correlation (Pearson)",
    annotation_col = annotation_col,
    annotation_colors = annotation_colors
)
dev.off()
log_message("INFO", paste("Correlation heatmap saved:", cor_heatmap_file))

# ------------------------------------------------------------------
# Sample Distance Heatmap
# ------------------------------------------------------------------
log_message("INFO", "Creating sample distance heatmap...")

# Calculate Euclidean distance between samples
sample_dist <- dist(t(log2_cpm), method = "euclidean")
sample_dist_matrix <- as.matrix(sample_dist)

# Generate distance heatmap
dist_heatmap_file <- file.path(plots_dir, "sample_distances.pdf")

pdf(dist_heatmap_file, width = 10, height = 8)
pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dist,
    clustering_distance_cols = sample_dist,
    clustering_method = "complete",
    color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
    main = "Sample-Sample Distance (Euclidean)",
    annotation_col = annotation_col,
    annotation_colors = annotation_colors
)
dev.off()
log_message("INFO", paste("Distance heatmap saved:", dist_heatmap_file))

# ==============================================================================
# COMPLETION
# ==============================================================================

log_message("INFO", "==========================================")
log_message("INFO", "Normalization completed successfully")
log_message("INFO", paste("Output directory:", opt$output))
log_message("INFO", "==========================================")
