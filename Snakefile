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


def _resolve_fastq_from_samplesheet(sample, read):
    """
    Resolve FASTQ path from explicit sample sheet columns when available.
    """
    # Prefer R1/R2 naming, but keep backward compatibility with forward/reverse.
    has_r = "R1" in samples_df.columns and "R2" in samples_df.columns
    has_forward_reverse = "forward" in samples_df.columns and "reverse" in samples_df.columns
    if not has_r and not has_forward_reverse:
        return None

    row = samples_df.loc[samples_df["sample"] == sample]
    if row.empty:
        return None

    if has_r:
        col = "R1" if str(read) == "1" else "R2"
    else:
        col = "forward" if str(read) == "1" else "reverse"
    value = row.iloc[0][col]
    if pd.isna(value) or str(value).strip() == "":
        return None

    return str(Path(config["data_dir"]) / str(value))


def _resolve_fastq_by_pattern(sample, read):
    """
    Resolve FASTQ path using common naming patterns, including Illumina _001 suffix.
    """
    data_dir = Path(config["data_dir"])
    read = str(read)
    patterns = [
        f"{sample}_R{read}.fastq.gz",
        f"{sample}_R{read}_001.fastq.gz",
        f"*{sample}*R{read}_001.fastq.gz",
        f"*{sample}*R{read}.fastq.gz",
    ]

    matches = []
    for pattern in patterns:
        matches.extend(sorted(data_dir.glob(pattern)))

    # Preserve order while removing duplicates.
    unique_matches = list(dict.fromkeys(matches))

    if len(unique_matches) == 1:
        return str(unique_matches[0])

    if len(unique_matches) == 0:
        raise ValueError(
            f"No FASTQ matched sample '{sample}' read R{read} in {data_dir}. "
            "Either set exact file names in samples.csv forward/reverse columns "
            "or rename files to a supported pattern."
        )

    raise ValueError(
        f"Multiple FASTQs matched sample '{sample}' read R{read}: "
        f"{[str(p) for p in unique_matches]}\n"
        "Please make sample identifiers more specific or set exact forward/reverse "
        "file names in samples.csv."
    )


def get_fastq_path(sample, read):
    """
    Resolve FASTQ path by preferring explicit sample-sheet file names, then patterns.
    """
    from_sheet = _resolve_fastq_from_samplesheet(sample, read)
    if from_sheet is not None:
        return from_sheet
    return _resolve_fastq_by_pattern(sample, read)

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
include: "workflow/rules/databases.smk"
include: "workflow/rules/preprocessing.smk" 
include: "workflow/rules/taxonomy.smk"
include: "workflow/rules/function.smk"
include: "workflow/rules/stats.smk"
include: "workflow/rules/assembly.smk"
include: "workflow/rules/binning.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/visualization.smk"