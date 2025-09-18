# Snakemake workflow for genomic data analysis pipeline
# 
# This file serves as the main entry point for Snakemake workflows.
# It integrates various genomic analysis workflows including DNA-seq,
# RNA-seq, single-cell RNA-seq, and ChIP-seq analysis.
#
# Author: Gabriel Demetrios Lafis
# Project: genomic-data-analysis-pipeline
# Created: 2025-09-18

# =============================================================================
# CONFIGURATION
# =============================================================================

# Configuration file
configfile: "config/config.yaml"

# =============================================================================
# TARGET RULES
# =============================================================================

# Default rule - entry point for the workflow
rule all:
    input:
        # Main analysis results
        "results/summary/analysis_complete.txt"

# =============================================================================
# WORKFLOW MODULES
# =============================================================================

# Include workflow modules
# Uncomment and modify paths as needed:
# include: "workflows/snakemake/dna_seq.smk"
# include: "workflows/snakemake/rna_seq.smk"
# include: "workflows/snakemake/scrna_seq.smk"
# include: "workflows/snakemake/chip_seq.smk"

# =============================================================================
# SUMMARY RULE
# =============================================================================

# Create summary file indicating workflow completion
rule create_summary:
    output:
        "results/summary/analysis_complete.txt"
    shell:
        """
        mkdir -p results/summary
        echo "Genomic analysis pipeline completed successfully" > {output}
        echo "Date: $(date)" >> {output}
        echo "Pipeline: genomic-data-analysis-pipeline" >> {output}
        """

# =============================================================================
# CLEAN RULES
# =============================================================================

# Clean temporary files
rule clean:
    shell:
        "rm -rf .snakemake logs results/temp"

# Clean all results
rule clean_all:
    shell:
        "rm -rf .snakemake logs results"
