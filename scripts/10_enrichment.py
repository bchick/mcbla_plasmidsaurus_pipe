#!/usr/bin/env python3
"""
10_enrichment.py - Perform functional enrichment analysis using GSEApy

Description:
    This script performs Gene Set Enrichment Analysis (GSEA) using GSEApy to
    identify biological pathways and processes enriched in differentially
    expressed genes. It uses pre-ranked GSEA with the MSigDB Hallmark gene set
    collection, which provides a curated set of gene sets representing well-defined
    biological states and processes.

Usage:
    python 10_enrichment.py --ranked <ranked_genes.rnk> --output <output_dir>
        [--organism human] [--gene-sets MSigDB_Hallmark_2020]

Arguments:
    --ranked        Path to ranked gene list (.rnk file from DE analysis) [required]
    --output        Output directory for enrichment results [required]
    --organism      Organism: human or mouse [default: human]
    --gene-sets     Gene set database(s), comma-separated [default: MSigDB_Hallmark_2020]
    --min-size      Minimum gene set size [default: 15]
    --max-size      Maximum gene set size [default: 500]
    --permutations  Number of permutations [default: 1000]
    --fdr           FDR threshold for significance [default: 0.25]
    --help          Display this help message

Dependencies:
    - Python >= 3.8
    - gseapy >= 0.12
    - pandas
    - matplotlib

Output:
    - gsea_report.txt               Summary of enriched gene sets
    - gsea_results_*.csv            Full results for each gene set database
    - plots/                        GSEA plots (enrichment plots, dot plots)

Example:
    python 10_enrichment.py \\
        --ranked results/09_de/ranked_genes_treatment_vs_control.rnk \\
        --output results/10_enrichment/ \\
        --organism human \\
        --gene-sets "MSigDB_Hallmark_2020,KEGG_2021_Human"

Author: Plasmidsaurus RNA-seq Pipeline
Date: 2024
Version: 1.0.0
"""

# ==============================================================================
# IMPORTS
# ==============================================================================

import os
import sys
import argparse
import logging
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

# ==============================================================================
# LOGGING SETUP
# ==============================================================================

def setup_logging() -> logging.Logger:
    """Configure logging with timestamp format."""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.INFO)

    formatter = logging.Formatter(
        '[%(asctime)s] [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger

logger = setup_logging()

# ==============================================================================
# CONSTANTS
# ==============================================================================

# Available gene set databases in Enrichr/GSEApy
# These are organism-specific where noted
AVAILABLE_GENE_SETS = {
    'human': [
        'MSigDB_Hallmark_2020',
        'KEGG_2021_Human',
        'GO_Biological_Process_2023',
        'GO_Molecular_Function_2023',
        'GO_Cellular_Component_2023',
        'Reactome_2022',
        'WikiPathway_2023_Human',
    ],
    'mouse': [
        'MSigDB_Hallmark_2020',
        'KEGG_2019_Mouse',
        'GO_Biological_Process_2023',
        'GO_Molecular_Function_2023',
        'GO_Cellular_Component_2023',
        'Reactome_2022',
        'WikiPathway_2023_Mouse',
    ]
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Perform functional enrichment analysis using GSEA',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --ranked ranked_genes.rnk --output results/enrichment/
  %(prog)s --ranked ranked_genes.rnk --output results/ --organism mouse
  %(prog)s --ranked ranked_genes.rnk --output results/ --gene-sets "KEGG_2021_Human,Reactome_2022"
        """
    )

    parser.add_argument(
        '--ranked', required=True,
        help='Path to ranked gene list (.rnk file)'
    )
    parser.add_argument(
        '--output', required=True,
        help='Output directory for enrichment results'
    )
    parser.add_argument(
        '--organism', default='human', choices=['human', 'mouse'],
        help='Organism for gene set selection (default: human)'
    )
    parser.add_argument(
        '--gene-sets', default='MSigDB_Hallmark_2020',
        help='Gene set database(s), comma-separated (default: MSigDB_Hallmark_2020)'
    )
    parser.add_argument(
        '--min-size', type=int, default=15,
        help='Minimum gene set size (default: 15)'
    )
    parser.add_argument(
        '--max-size', type=int, default=500,
        help='Maximum gene set size (default: 500)'
    )
    parser.add_argument(
        '--permutations', type=int, default=1000,
        help='Number of permutations for GSEA (default: 1000)'
    )
    parser.add_argument(
        '--fdr', type=float, default=0.25,
        help='FDR threshold for significance (default: 0.25)'
    )

    return parser.parse_args()


def load_ranked_genes(filepath: str) -> pd.Series:
    """
    Load ranked gene list from .rnk file.

    The .rnk format is a two-column tab-separated file:
    - Column 1: Gene identifier (symbol or ID)
    - Column 2: Ranking metric (e.g., signed -log10 p-value)

    Args:
        filepath: Path to .rnk file

    Returns:
        Pandas Series with gene names as index and rank metric as values
    """
    logger.info(f"Loading ranked gene list: {filepath}")

    if not os.path.exists(filepath):
        logger.error(f"Ranked gene file not found: {filepath}")
        sys.exit(1)

    # Read the ranked gene list
    # .rnk files have no header
    df = pd.read_csv(filepath, sep='\t', header=None, names=['gene', 'rank'])

    # Remove any duplicate genes (keep first occurrence)
    if df['gene'].duplicated().any():
        n_dups = df['gene'].duplicated().sum()
        logger.warning(f"Removing {n_dups} duplicate gene entries")
        df = df.drop_duplicates(subset='gene', keep='first')

    # Convert to Series
    ranked_genes = df.set_index('gene')['rank']

    # Sort by rank (descending - most upregulated first)
    ranked_genes = ranked_genes.sort_values(ascending=False)

    logger.info(f"Loaded {len(ranked_genes)} genes")
    logger.info(f"Rank range: {ranked_genes.min():.2f} to {ranked_genes.max():.2f}")

    return ranked_genes


def run_gsea_prerank(
    ranked_genes: pd.Series,
    gene_set: str,
    output_dir: str,
    min_size: int = 15,
    max_size: int = 500,
    permutations: int = 1000
) -> Optional[pd.DataFrame]:
    """
    Run pre-ranked GSEA analysis.

    Pre-ranked GSEA uses a pre-computed ranking metric (e.g., from differential
    expression analysis) rather than computing rankings from expression data.
    This approach is useful when:
    1. The experimental design is complex
    2. You want to use a custom ranking metric
    3. You want to analyze publicly available ranked gene lists

    Args:
        ranked_genes: Series with gene names as index and rank metric as values
        gene_set: Name of gene set database to use
        output_dir: Directory for GSEA output
        min_size: Minimum gene set size to consider
        max_size: Maximum gene set size to consider
        permutations: Number of permutations for p-value estimation

    Returns:
        DataFrame with GSEA results, or None if analysis failed
    """
    logger.info(f"Running GSEA with gene set: {gene_set}")

    try:
        # Run pre-ranked GSEA
        # GSEApy's prerank function performs the core GSEA algorithm:
        # 1. Walks down the ranked list
        # 2. Increases running sum when gene is in set, decreases otherwise
        # 3. Enrichment score (ES) is maximum deviation from zero
        # 4. Normalized ES (NES) accounts for gene set size
        # 5. P-value estimated by permutation testing

        gsea_results = gp.prerank(
            rnk=ranked_genes,
            gene_sets=gene_set,
            outdir=os.path.join(output_dir, gene_set.replace(' ', '_')),
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutations,
            # Use 'signal_to_noise' weighting for ranking metric
            weighted_score_type=1,
            # Number of processes (set to 1 for reproducibility)
            processes=4,
            # Seed for reproducibility
            seed=42,
            # Verbose output
            verbose=True
        )

        # Extract results DataFrame
        results_df = gsea_results.res2d

        if results_df is not None and len(results_df) > 0:
            logger.info(f"GSEA completed: {len(results_df)} gene sets tested")
            return results_df
        else:
            logger.warning(f"No results returned for gene set: {gene_set}")
            return None

    except Exception as e:
        logger.error(f"GSEA failed for {gene_set}: {str(e)}")
        return None


def generate_summary_report(
    all_results: dict,
    output_dir: str,
    fdr_threshold: float = 0.25
) -> None:
    """
    Generate summary report of enrichment results.

    Args:
        all_results: Dictionary mapping gene set names to result DataFrames
        output_dir: Output directory
        fdr_threshold: FDR threshold for reporting significant results
    """
    logger.info("Generating summary report...")

    report_path = os.path.join(output_dir, "gsea_report.txt")

    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("GENE SET ENRICHMENT ANALYSIS REPORT\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("=" * 80 + "\n\n")

        for gene_set_name, results_df in all_results.items():
            if results_df is None:
                f.write(f"\n{gene_set_name}: No results\n")
                continue

            f.write(f"\n{'=' * 80}\n")
            f.write(f"Gene Set Database: {gene_set_name}\n")
            f.write(f"{'=' * 80}\n\n")

            # Filter significant results
            sig_up = results_df[
                (results_df['FDR q-val'] < fdr_threshold) &
                (results_df['NES'] > 0)
            ].sort_values('NES', ascending=False)

            sig_down = results_df[
                (results_df['FDR q-val'] < fdr_threshold) &
                (results_df['NES'] < 0)
            ].sort_values('NES', ascending=True)

            f.write(f"Total gene sets tested: {len(results_df)}\n")
            f.write(f"Significant (FDR < {fdr_threshold}): {len(sig_up) + len(sig_down)}\n")
            f.write(f"  - Upregulated: {len(sig_up)}\n")
            f.write(f"  - Downregulated: {len(sig_down)}\n\n")

            if len(sig_up) > 0:
                f.write("TOP UPREGULATED PATHWAYS:\n")
                f.write("-" * 60 + "\n")
                for idx, row in sig_up.head(10).iterrows():
                    f.write(f"  {row['Term']}\n")
                    f.write(f"    NES: {row['NES']:.3f}, FDR: {row['FDR q-val']:.4f}\n")
                f.write("\n")

            if len(sig_down) > 0:
                f.write("TOP DOWNREGULATED PATHWAYS:\n")
                f.write("-" * 60 + "\n")
                for idx, row in sig_down.head(10).iterrows():
                    f.write(f"  {row['Term']}\n")
                    f.write(f"    NES: {row['NES']:.3f}, FDR: {row['FDR q-val']:.4f}\n")
                f.write("\n")

    logger.info(f"Summary report saved: {report_path}")


def create_dot_plot(
    results_df: pd.DataFrame,
    gene_set_name: str,
    output_dir: str,
    top_n: int = 20,
    fdr_threshold: float = 0.25
) -> None:
    """
    Create dot plot visualization of enrichment results.

    Args:
        results_df: GSEA results DataFrame
        gene_set_name: Name of gene set database
        output_dir: Output directory for plots
        top_n: Number of top pathways to show
        fdr_threshold: FDR threshold for coloring
    """
    if results_df is None or len(results_df) == 0:
        return

    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Select top pathways by absolute NES
    top_results = results_df.nlargest(top_n, 'NES', keep='first')
    bottom_results = results_df.nsmallest(top_n, 'NES', keep='first')
    plot_data = pd.concat([top_results, bottom_results]).drop_duplicates()

    if len(plot_data) == 0:
        return

    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(8, len(plot_data) * 0.3)))

    # Create dot plot
    scatter = ax.scatter(
        plot_data['NES'],
        range(len(plot_data)),
        c=-np.log10(plot_data['FDR q-val'].clip(1e-10, 1)),
        s=plot_data['Tag %'].str.rstrip('%').astype(float) * 3 if 'Tag %' in plot_data.columns else 100,
        cmap='RdYlBu_r',
        alpha=0.7
    )

    # Customize plot
    ax.set_yticks(range(len(plot_data)))
    ax.set_yticklabels(plot_data['Term'].values, fontsize=8)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=10)
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, label='-log10(FDR)')

    # Title
    ax.set_title(f'GSEA Results: {gene_set_name}', fontsize=12, fontweight='bold')

    plt.tight_layout()

    # Save plot
    plot_path = os.path.join(plots_dir, f"dotplot_{gene_set_name.replace(' ', '_')}.pdf")
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Dot plot saved: {plot_path}")


# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

def main():
    """Main function to run enrichment analysis."""
    logger.info("=" * 50)
    logger.info("Starting Functional Enrichment Analysis")
    logger.info("=" * 50)

    # Parse arguments
    args = parse_arguments()

    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    os.makedirs(os.path.join(args.output, "plots"), exist_ok=True)

    # Log parameters
    logger.info(f"Ranked gene file: {args.ranked}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Organism: {args.organism}")
    logger.info(f"Gene sets: {args.gene_sets}")
    logger.info(f"Gene set size: {args.min_size}-{args.max_size}")
    logger.info(f"Permutations: {args.permutations}")
    logger.info(f"FDR threshold: {args.fdr}")

    # Load ranked genes
    ranked_genes = load_ranked_genes(args.ranked)

    # Parse gene sets to analyze
    gene_sets = [gs.strip() for gs in args.gene_sets.split(',')]

    # Run GSEA for each gene set database
    all_results = {}

    for gene_set in gene_sets:
        logger.info(f"\n{'=' * 50}")
        logger.info(f"Analyzing gene set: {gene_set}")
        logger.info("=" * 50)

        results = run_gsea_prerank(
            ranked_genes=ranked_genes,
            gene_set=gene_set,
            output_dir=args.output,
            min_size=args.min_size,
            max_size=args.max_size,
            permutations=args.permutations
        )

        all_results[gene_set] = results

        if results is not None:
            # Save full results
            results_path = os.path.join(
                args.output,
                f"gsea_results_{gene_set.replace(' ', '_')}.csv"
            )
            results.to_csv(results_path, index=False)
            logger.info(f"Results saved: {results_path}")

            # Count significant results
            n_sig = len(results[results['FDR q-val'] < args.fdr])
            logger.info(f"Significant pathways (FDR < {args.fdr}): {n_sig}")

            # Create dot plot
            try:
                import numpy as np
                create_dot_plot(
                    results,
                    gene_set,
                    args.output,
                    fdr_threshold=args.fdr
                )
            except Exception as e:
                logger.warning(f"Could not create dot plot: {e}")

    # Generate summary report
    generate_summary_report(all_results, args.output, args.fdr)

    logger.info("\n" + "=" * 50)
    logger.info("Functional enrichment analysis completed")
    logger.info(f"Output directory: {args.output}")
    logger.info("=" * 50)


if __name__ == "__main__":
    main()
