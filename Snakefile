"""
Metagenomics Shotgun Analysis Pipeline
======================================
A minimal Snakemake pipeline for metagenomics shotgun sequencing analysis

Tools used:
- FastQC: Quality control
- KneadData: Host removal and quality filtering
- MetaPhlAn4: Taxonomic profiling
- HUMAnN3: Functional profiling
- MaAsLin2: Statistical analysis
"""

import pandas as pd
from pathlib import Path

# Load configuration
configfile: "config/config.yaml"

# Load sample sheet
samples_df = pd.read_csv(config["sample_sheet"])
SAMPLES = samples_df["sample"].tolist()

# Define all output files
rule all:
    input:
        # QC reports
        expand("results/fastqc/{sample}_R{read}_fastqc.html", sample=SAMPLES, read=[1,2]),
        # Host removal
        expand("results/kneaddata/{sample}_kneaddata_paired_1.fastq", sample=SAMPLES),
        expand("results/kneaddata/{sample}_kneaddata_paired_2.fastq", sample=SAMPLES),
        # Taxonomic profiling
        expand("results/metaphlan/{sample}_metaphlan.txt", sample=SAMPLES),
        "results/metaphlan/merged_abundance_table.txt",
        # Functional profiling
        expand("results/humann/{sample}_genefamilies.tsv", sample=SAMPLES),
        expand("results/humann/{sample}_pathabundance.tsv", sample=SAMPLES),
        expand("results/humann/{sample}_pathcoverage.tsv", sample=SAMPLES),
        "results/humann/merged_genefamilies.tsv",
        "results/humann/merged_pathabundance.tsv",
        # Basic visualization
        "results/visualization/pca_taxonomy.png",
        "results/visualization/taxonomy_barplot.png",
        # Final report
        "results/multiqc_report.html"

# Advanced analysis (optional - run with specific targets)
rule all_advanced:
    input:
        # Assembly
        expand("results/assembly/{sample}/final.contigs.fa", sample=SAMPLES),
        # Binning
        expand("results/binning/{sample}/checkm/quality_summary.tsv", sample=SAMPLES),
        # Annotation 
        expand("results/annotation/{sample}/dram/genome_summaries.tsv", sample=SAMPLES),
        # Advanced visualization
        "results/visualization/metagenomics_dashboard.html"

# Include workflow rules
include: "workflow/rules/qc.smk"
include: "workflow/rules/preprocessing.smk" 
include: "workflow/rules/taxonomy.smk"
include: "workflow/rules/function.smk"
include: "workflow/rules/stats.smk"
include: "workflow/rules/assembly.smk"
include: "workflow/rules/binning.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/visualization.smk"