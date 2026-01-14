#!/usr/bin/env Rscript
#
# 09_differential_expression.R - Perform differential expression analysis using edgeR
#
# Description:
#   This script performs differential expression (DE) analysis using edgeR's
#   quasi-likelihood framework. It tests for genes with significant expression
#   changes between conditions while controlling for false discovery rate (FDR).
#   The script supports multiple contrasts and generates comprehensive output
#   including MA plots, volcano plots, and ranked gene lists for GSEA.
#
# Usage:
#   Rscript 09_differential_expression.R --dge <dge_object.rds> --metadata <sample_info> \
#       --output <output_dir> --contrast <treatment-control>
#
# Arguments:
#   --dge         Path to normalized DGEList RDS file from 08_normalize.R (required)
#   --metadata    Path to sample metadata TSV file (required)
#   --output      Output directory for DE results (required)
#   --contrast    Contrast to test, format: "condition1-condition2" (required)
#   --fdr         FDR threshold for significance (default: 0.05)
#   --lfc         Log2 fold change threshold (default: 1.0)
#   --help        Display this help message
#
# Dependencies:
#   - R >= 4.0
#   - edgeR >= 4.0.16
#   - ggplot2
#   - ggrepel
#
# Output:
#   - de_results_{contrast}.txt       Full DE results table
#   - de_significant_{contrast}.txt   Significant DE genes only
#   - ranked_genes_{contrast}.rnk     Ranked gene list for GSEA
#   - plots/ma_plot_{contrast}.pdf    MA plot
#   - plots/volcano_{contrast}.pdf    Volcano plot
#
# Example:
#   Rscript 09_differential_expression.R \
#       --dge results/08_normalization/dge_normalized.rds \
#       --metadata config/samples.tsv \
#       --output results/09_de/ \
#       --contrast "treatment-control"
#
# Author: Plasmidsaurus RNA-seq Pipeline
# Date: 2024
# Version: 1.0.0

# ==============================================================================
# SETUP AND LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(edgeR)
    library(ggplot2)
    library(ggrepel)
})

set.seed(42)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

log_message <- function(level, message) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", timestamp, level, message), file = stderr())
}

die <- function(message) {
    log_message("ERROR", message)
    quit(status = 1)
}

#' Parse contrast string into condition names
#'
#' @param contrast_str Contrast string in format "condition1-condition2"
#' @return Named list with 'treatment' and 'control' conditions
parse_contrast <- function(contrast_str) {
    parts <- strsplit(contrast_str, "-")[[1]]
    if (length(parts) != 2) {
        die(paste("Invalid contrast format:", contrast_str,
                  "- Expected format: condition1-condition2"))
    }
    list(treatment = trimws(parts[1]), control = trimws(parts[2]))
}

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

option_list <- list(
    make_option(c("--dge"), type = "character", default = NULL,
                help = "Path to normalized DGEList RDS file [required]"),
    make_option(c("--metadata"), type = "character", default = NULL,
                help = "Path to sample metadata TSV file [required]"),
    make_option(c("--output"), type = "character", default = NULL,
                help = "Output directory for DE results [required]"),
    make_option(c("--contrast"), type = "character", default = NULL,
                help = "Contrast to test: condition1-condition2 [required]"),
    make_option(c("--fdr"), type = "double", default = 0.05,
                help = "FDR threshold for significance [default: %default]"),
    make_option(c("--lfc"), type = "double", default = 1.0,
                help = "Log2 fold change threshold [default: %default]")
)

opt_parser <- OptionParser(
    option_list = option_list,
    description = "Perform differential expression analysis using edgeR"
)
opt <- parse_args(opt_parser)

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

log_message("INFO", "==========================================")
log_message("INFO", "Starting Differential Expression Analysis")
log_message("INFO", "==========================================")

if (is.null(opt$dge)) die("DGEList RDS file is required (--dge)")
if (is.null(opt$metadata)) die("Sample metadata file is required (--metadata)")
if (is.null(opt$output)) die("Output directory is required (--output)")
if (is.null(opt$contrast)) die("Contrast is required (--contrast)")

if (!file.exists(opt$dge)) die(paste("DGEList file not found:", opt$dge))
if (!file.exists(opt$metadata)) die(paste("Metadata file not found:", opt$metadata))

dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(opt$output, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Parse contrast
contrast <- parse_contrast(opt$contrast)
contrast_name <- gsub("-", "_vs_", opt$contrast)

log_message("INFO", paste("Contrast:", contrast$treatment, "vs", contrast$control))
log_message("INFO", paste("FDR threshold:", opt$fdr))
log_message("INFO", paste("LFC threshold:", opt$lfc))

# ==============================================================================
# LOAD DATA
# ==============================================================================

log_message("INFO", "Loading normalized DGEList...")
dge <- readRDS(opt$dge)

log_message("INFO", "Loading sample metadata...")
metadata <- read.table(opt$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Verify condition column exists
if (!"condition" %in% colnames(metadata)) {
    die("Metadata must contain 'condition' column")
}

# Verify contrast conditions exist in metadata
available_conditions <- unique(metadata$condition)
if (!contrast$treatment %in% available_conditions) {
    die(paste("Treatment condition not found:", contrast$treatment,
              "- Available:", paste(available_conditions, collapse = ", ")))
}
if (!contrast$control %in% available_conditions) {
    die(paste("Control condition not found:", contrast$control,
              "- Available:", paste(available_conditions, collapse = ", ")))
}

# Match metadata to DGE sample order
metadata <- metadata[match(colnames(dge), metadata$sample_id), ]

log_message("INFO", paste("Samples in analysis:", ncol(dge)))
log_message("INFO", paste("Genes in analysis:", nrow(dge)))

# ==============================================================================
# DESIGN MATRIX
# ==============================================================================

log_message("INFO", "Creating design matrix...")

# Create factor for condition with control as reference level
# This ensures coefficients represent fold change relative to control
condition <- factor(metadata$condition)
condition <- relevel(condition, ref = contrast$control)

# Create design matrix
# ~ 0 + condition creates a cell means model where each coefficient is a group mean
# This makes it easier to specify arbitrary contrasts
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

log_message("INFO", "Design matrix:")
print(design)

# ==============================================================================
# DISPERSION ESTIMATION
# ==============================================================================

log_message("INFO", "Estimating dispersion parameters...")

# edgeR models gene counts using the Negative Binomial distribution
# which has two parameters: mean and dispersion
# Dispersion captures biological variability (how much counts vary between replicates)

# estimateDisp performs:
# 1. Common dispersion: single dispersion value for all genes
# 2. Trended dispersion: dispersion as function of abundance (CPM)
# 3. Tagwise dispersion: gene-specific dispersion, shrunk toward trend

dge <- estimateDisp(dge, design, robust = TRUE)

log_message("INFO", sprintf("Common dispersion: %.4f", dge$common.dispersion))
log_message("INFO", sprintf("Common BCV (biological coefficient of variation): %.2f%%",
                            sqrt(dge$common.dispersion) * 100))

# ==============================================================================
# QUASI-LIKELIHOOD MODEL FITTING
# ==============================================================================

log_message("INFO", "Fitting quasi-likelihood model...")

# The quasi-likelihood (QL) framework provides more robust inference than
# standard likelihood methods, especially with small sample sizes
#
# glmQLFit fits a GLM with quasi-likelihood dispersion estimation
# The QL approach:
# 1. Adds an extra layer of variability modeling (QL dispersion)
# 2. Accounts for uncertainty in dispersion estimation
# 3. Provides more accurate control of false discovery rates

fit <- glmQLFit(dge, design, robust = TRUE)

# ==============================================================================
# DIFFERENTIAL EXPRESSION TESTING
# ==============================================================================

log_message("INFO", paste("Testing contrast:", opt$contrast))

# Create contrast vector
# This specifies which coefficients to compare
# makeContrasts creates the vector automatically from a formula
contrast_vector <- makeContrasts(
    contrasts = paste0(contrast$treatment, "-", contrast$control),
    levels = design
)

# Perform quasi-likelihood F-test
# The QL F-test is preferred over the LRT (likelihood ratio test) because:
# 1. It accounts for uncertainty in dispersion estimation
# 2. It provides better FDR control with small sample sizes
# 3. It's more robust to outliers

qlf_result <- glmQLFTest(fit, contrast = contrast_vector)

# Extract results
# topTags returns a table with logFC, logCPM, F-statistic, p-value, and FDR
de_results <- topTags(qlf_result, n = Inf, sort.by = "PValue")$table

# Add gene information if available
if (!is.null(dge$genes)) {
    de_results <- cbind(
        gene_id = rownames(de_results),
        dge$genes[rownames(de_results), , drop = FALSE],
        de_results
    )
} else {
    de_results <- cbind(gene_id = rownames(de_results), de_results)
}

# ==============================================================================
# SUMMARIZE RESULTS
# ==============================================================================

log_message("INFO", "Summarizing differential expression results...")

# Apply significance thresholds
is_significant <- de_results$FDR < opt$fdr & abs(de_results$logFC) >= opt$lfc

de_results$significant <- ifelse(is_significant, "Yes", "No")
de_results$direction <- ifelse(de_results$logFC > 0, "Up", "Down")
de_results$direction[!is_significant] <- "NS"

# Count significant genes
n_up <- sum(de_results$direction == "Up", na.rm = TRUE)
n_down <- sum(de_results$direction == "Down", na.rm = TRUE)
n_total <- n_up + n_down

log_message("INFO", "==========================================")
log_message("INFO", "Differential Expression Summary")
log_message("INFO", "==========================================")
log_message("INFO", sprintf("Total genes tested: %d", nrow(de_results)))
log_message("INFO", sprintf("Significant genes (FDR < %.2f, |LFC| >= %.1f): %d",
                            opt$fdr, opt$lfc, n_total))
log_message("INFO", sprintf("  Upregulated in %s: %d", contrast$treatment, n_up))
log_message("INFO", sprintf("  Downregulated in %s: %d", contrast$treatment, n_down))

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

log_message("INFO", "Saving results...")

# Save full results
full_output <- file.path(opt$output, paste0("de_results_", contrast_name, ".txt"))
write.table(de_results, full_output, sep = "\t", quote = FALSE, row.names = FALSE)
log_message("INFO", paste("Full results saved:", full_output))

# Save significant genes only
sig_results <- de_results[de_results$significant == "Yes", ]
sig_output <- file.path(opt$output, paste0("de_significant_", contrast_name, ".txt"))
write.table(sig_results, sig_output, sep = "\t", quote = FALSE, row.names = FALSE)
log_message("INFO", paste("Significant genes saved:", sig_output))

# Create ranked gene list for GSEA
# GSEA uses pre-ranked gene lists based on a ranking metric
# We use signed -log10(p-value) which captures both significance and direction
rank_metric <- sign(de_results$logFC) * -log10(de_results$PValue)
ranked_genes <- data.frame(
    gene = de_results$gene_id,
    rank = rank_metric
)
ranked_genes <- ranked_genes[order(ranked_genes$rank, decreasing = TRUE), ]

rnk_output <- file.path(opt$output, paste0("ranked_genes_", contrast_name, ".rnk"))
write.table(ranked_genes, rnk_output, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
log_message("INFO", paste("Ranked gene list saved:", rnk_output))

# ==============================================================================
# GENERATE PLOTS
# ==============================================================================

log_message("INFO", "Generating visualization plots...")

# ------------------------------------------------------------------
# MA Plot
# ------------------------------------------------------------------
log_message("INFO", "Creating MA plot...")

# MA plot shows log fold change vs average expression
# Helps identify expression-dependent bias and visualize DE pattern
ma_df <- data.frame(
    logCPM = de_results$logCPM,
    logFC = de_results$logFC,
    significant = factor(de_results$direction, levels = c("Up", "Down", "NS"))
)

ma_plot <- ggplot(ma_df, aes(x = logCPM, y = logFC, color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = c(-opt$lfc, 0, opt$lfc),
               linetype = c("dashed", "solid", "dashed"),
               color = c("blue", "black", "red"), alpha = 0.5) +
    scale_color_manual(
        values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "gray60"),
        labels = c(
            "Up" = sprintf("Up (%d)", n_up),
            "Down" = sprintf("Down (%d)", n_down),
            "NS" = "Not significant"
        )
    ) +
    labs(
        title = paste("MA Plot:", contrast$treatment, "vs", contrast$control),
        x = "Average log2 CPM",
        y = "log2 Fold Change",
        color = "Regulation"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom"
    )

ma_file <- file.path(plots_dir, paste0("ma_plot_", contrast_name, ".pdf"))
ggsave(ma_file, ma_plot, width = 10, height = 8)
log_message("INFO", paste("MA plot saved:", ma_file))

# ------------------------------------------------------------------
# Volcano Plot
# ------------------------------------------------------------------
log_message("INFO", "Creating volcano plot...")

# Volcano plot shows -log10(p-value) vs log fold change
# Creates a "volcano" shape with significant genes at the top corners
volcano_df <- data.frame(
    logFC = de_results$logFC,
    negLog10P = -log10(de_results$PValue),
    significant = factor(de_results$direction, levels = c("Up", "Down", "NS")),
    gene_id = de_results$gene_id
)

# Identify top genes to label
top_genes <- head(de_results[order(de_results$PValue), ], 20)
volcano_df$label <- ifelse(volcano_df$gene_id %in% top_genes$gene_id,
                           as.character(volcano_df$gene_id), "")

# Use gene_name if available
if ("gene_name" %in% colnames(de_results)) {
    gene_name_map <- setNames(de_results$gene_name, de_results$gene_id)
    volcano_df$label <- ifelse(volcano_df$label != "",
                               gene_name_map[volcano_df$label],
                               "")
}

volcano_plot <- ggplot(volcano_df, aes(x = logFC, y = negLog10P, color = significant)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_vline(xintercept = c(-opt$lfc, opt$lfc), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(opt$fdr), linetype = "dashed", alpha = 0.5) +
    geom_text_repel(
        aes(label = label),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5
    ) +
    scale_color_manual(
        values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "gray60"),
        labels = c(
            "Up" = sprintf("Up (%d)", n_up),
            "Down" = sprintf("Down (%d)", n_down),
            "NS" = "Not significant"
        )
    ) +
    labs(
        title = paste("Volcano Plot:", contrast$treatment, "vs", contrast$control),
        x = "log2 Fold Change",
        y = "-log10(P-value)",
        color = "Regulation"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom"
    )

volcano_file <- file.path(plots_dir, paste0("volcano_", contrast_name, ".pdf"))
ggsave(volcano_file, volcano_plot, width = 10, height = 8)
log_message("INFO", paste("Volcano plot saved:", volcano_file))

# ==============================================================================
# COMPLETION
# ==============================================================================

log_message("INFO", "==========================================")
log_message("INFO", "Differential expression analysis completed")
log_message("INFO", paste("Output directory:", opt$output))
log_message("INFO", "==========================================")
