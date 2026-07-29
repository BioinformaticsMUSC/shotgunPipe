# Metagenomics Shotgun Analysis Pipeline

A minimal Snakemake pipeline for metagenomics shotgun sequencing analysis, implementing best practices for microbiome data processing.

## Pipeline Overview

This pipeline processes paired-end shotgun metagenomic sequencing data through the following steps:

**Core Analysis:**
1. **Quality Control** (FastQC)
2. **Host Removal** (KneadData)  
3. **Taxonomic Profiling** (MetaPhlAn4)
4. **Functional Profiling** (HUMAnN3)
5. **Statistical Analysis** (MaAsLin2)
6. **Quality Report** (MultiQC)

**Advanced Analysis (Optional):**
7. **Assembly** (MEGAHIT/metaSPAdes)
8. **Binning** (MetaBAT2/CONCOCT/DAS Tool)
9. **MAG Annotation** (Prokka/DRAM/eggNOG-mapper)
10. **Visualization** (Krona/PCA/Diversity plots)

## Tools Used

**Core Tools:**
- **FastQC**: Quality control assessment of raw sequencing data
- **KneadData**: Host contamination removal and quality filtering
- **MetaPhlAn4**: Species-level taxonomic profiling
- **HUMAnN3**: Functional profiling and pathway analysis
- **MaAsLin2**: Multivariable statistical analysis for microbiome data
- **MultiQC**: Aggregate reporting of QC metrics

**Advanced Tools (Optional):**
- **MEGAHIT/metaSPAdes**: De novo metagenome assembly
- **MetaBAT2/CONCOCT**: Genome binning for MAG recovery
- **DAS Tool**: Bin optimization and dereplication
- **CheckM**: MAG quality assessment
- **Prokka**: Rapid prokaryotic genome annotation
- **DRAM**: Comprehensive metabolic annotation
- **eggNOG-mapper**: Functional annotation using eggNOG database
- **Krona**: Interactive taxonomic visualization
- **Custom scripts**: PCA, diversity analysis, and visualization

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- [Snakemake](https://snakemake.readthedocs.io/)

### Install Snakemake

```bash
# Install snakemake
conda install -c bioconda -c conda-forge snakemake
```

### Download Reference Databases

You'll need to download and configure the following databases:

The easiest path is to let Snakemake prepare them all:

```bash
snakemake setup_databases --use-conda --cores 1
```

This fills `resources/databases/` with the host genome, MetaPhlAn, HUMAnN, eggNOG, and dbCAN databases used by the workflow.

#### 1. Host Genome (e.g., Human GRCh38)
```bash
# Download human genome for KneadData
kneaddata_database --download human_genome bowtie2 /path/to/host/genome/
```

#### 2. MetaPhlAn4 Database
```bash
# Download MetaPhlAn database
metaphlan --install --bowtie2db /path/to/metaphlan/db/
```

#### 3. HUMAnN3 Databases
```bash
# Download ChocoPhlAn database
humann_databases --download chocophlan full /path/to/chocophlan/

# Download UniRef90 database
humann_databases --download uniref uniref90_diamond /path/to/uniref90/
```

## Configuration

### 1. Update Database Paths

Edit `config/config.yaml` to specify the correct paths to your reference databases:

```yaml
db_root: "resources/databases"

databases:
  host_genome: "resources/databases/kneaddata/human"
  metaphlan_db: "resources/databases/metaphlan"
  chocophlan_db: "resources/databases/humann/chocophlan"
  uniref_db: "resources/databases/humann/uniref90"

eggnog:
  database_dir: "resources/databases/eggnog"

cazy:
  database_dir: "resources/databases/dbcan"
```

### 2. Prepare Sample Sheet

Edit `config/samples.csv` with your sample information:

```csv
sample,sample_id,forward,reverse,condition
sample01,S001,sample01_R1.fastq.gz,sample01_R2.fastq.gz,treatment
sample02,S002,sample02_R1.fastq.gz,sample02_R2.fastq.gz,control
```

### 3. Data Organization

Organize your raw data as follows:

```
shotgun_pipe/
├── data/
│   └── raw/
│       ├── sample01_R1.fastq.gz
│       ├── sample01_R2.fastq.gz
│       ├── sample02_R1.fastq.gz
│       └── sample02_R2.fastq.gz
```

## Usage

### Run Complete Pipeline

```bash
# Dry run to check workflow
snakemake --dry-run --cores 1

# Run core pipeline (basic analysis)
snakemake --use-conda --cores 8

# Run advanced pipeline (includes assembly, binning, annotation)
snakemake all_advanced --use-conda --cores 8

# Run with specific number of jobs
snakemake --use-conda --jobs 4
```

### Run Specific Modules

```bash
# Core analysis only
snakemake results/multiqc_report.html --use-conda --cores 4

# Taxonomic profiling only
snakemake results/metaphlan/merged_abundance_table.txt --use-conda --cores 4

# Assembly only
snakemake results/assembly/all_contigs_combined.fa --use-conda --cores 8

# Binning only (requires assembly)
snakemake results/binning/sample01/checkm/quality_summary.tsv --use-conda --cores 8

# Visualization only
snakemake results/visualization/metagenomics_dashboard.html --use-conda --cores 4
```

### Cluster Execution

For SLURM clusters:

```bash
snakemake --use-conda --cluster "sbatch -t 24:00:00 -c {threads} --mem={resources.mem_mb}" --jobs 10
```

## Output Structure

```
results/
├── fastqc/                    # FastQC reports
├── kneaddata/                 # Host-removed reads
├── metaphlan/                 # Taxonomic profiles
│   └── merged_abundance_table.txt
├── humann/                    # Functional profiles
│   ├── merged_genefamilies_relab.tsv
│   └── merged_pathabundance_relab.tsv
├── maaslin2/                  # Statistical analysis
│   ├── taxonomy/
│   └── function/
├── assembly/                  # Assembled contigs (optional)
│   ├── sample01/final.contigs.fa
│   └── all_contigs_combined.fa
├── binning/                   # MAG binning results (optional)
│   ├── sample01/
│   │   ├── bins/             # MetaBAT2 bins
│   │   ├── concoct_bins/     # CONCOCT bins
│   │   ├── dastool_bins/     # Optimized bins
│   │   └── checkm/           # Quality assessment
├── annotation/                # MAG annotations (optional)
│   ├── sample01/
│   │   ├── prokka/           # Prokka annotations
│   │   ├── dram/             # DRAM annotations
│   │   └── eggnog/           # eggNOG annotations
├── visualization/             # Plots and visualizations
│   ├── pca_taxonomy.png
│   ├── taxonomy_barplot.png
│   ├── metagenomics_dashboard.html
│   └── *_krona.html
└── multiqc_report.html        # Summary report
```

## Key Output Files

**Core Analysis:**
- **`multiqc_report.html`**: Comprehensive QC summary
- **`metaphlan/merged_abundance_table.txt`**: Species abundance table
- **`humann/merged_pathabundance_relab.tsv`**: Pathway abundance table
- **`humann/merged_genefamilies_relab.tsv`**: Gene family abundance table
- **`maaslin2/taxonomy/`**: Statistical results for taxonomic data
- **`maaslin2/function/`**: Statistical results for functional data
- **`visualization/pca_taxonomy.png`**: PCA plot of samples
- **`visualization/taxonomy_barplot.png`**: Taxonomic composition plot

**Advanced Analysis:**
- **`assembly/all_contigs_combined.fa`**: Assembled contigs from all samples
- **`binning/{sample}/checkm/quality_summary.tsv`**: MAG quality metrics
- **`annotation/{sample}/dram/genome_summaries.tsv`**: Comprehensive MAG annotations
- **`visualization/metagenomics_dashboard.html`**: Interactive analysis dashboard
- **`visualization/*_krona.html`**: Interactive taxonomic trees

## Customization

### Adjust Tool Parameters

Modify `config/config.yaml` to adjust tool-specific parameters:

```yaml
kneaddata:
  threads: 8
  trimmomatic_options: "SLIDINGWINDOW:4:20 MINLEN:50"
  bowtie2_options: "--very-sensitive --dovetail"

metaphlan:
  threads: 4
  add_viruses: true
  unknown_estimation: false

humann:
  threads: 8
  search_mode: "uniref90"
  bypass_translated_search: false

maaslin2:
  min_abundance: 0.0
  min_prevalence: 0.1
  normalization: "TSS"
  transform: "LOG"
  analysis_method: "LM"
```

### Add Custom Rules

Create additional rules in the `workflow/rules/` directory and include them in the main `Snakefile`.

## Troubleshooting

### Common Issues

1. **Database paths**: Ensure all database paths in `config.yaml` are correct
2. **Sample names**: Check that sample names in `samples.csv` match your fastq file names
3. **Memory requirements**: Increase memory allocation for resource-intensive steps
4. **Environment conflicts**: Use `--use-conda` to ensure proper environment isolation

### Resource Requirements

**Core Analysis:**
- **RAM**: 16-32 GB recommended
- **Disk**: ~50-100 GB for intermediate files
- **CPU**: 8-16 cores for optimal performance

**Advanced Analysis (with assembly/binning):**
- **RAM**: 64-128 GB recommended (especially for metaSPAdes)
- **Disk**: ~200-500 GB for intermediate files (depending on dataset size)
- **CPU**: 16-32 cores for optimal performance

## Citation

If you use this pipeline, please cite the following tools:

- **KneadData**: https://github.com/biobakery/kneaddata
- **MetaPhlAn4**: Blanco-Míguez et al. eLife 2023
- **HUMAnN3**: Franzosa et al. Nature Methods 2018
- **MaAsLin2**: Mallick et al. Nature Methods 2021
- **Snakemake**: Köster & Rahmann, Bioinformatics 2012

## Support

For pipeline-specific issues, please check:
1. Snakemake logs in `.snakemake/log/`
2. Tool-specific logs in respective output directories
3. Environment setup with `conda list` in each environment

## License

This pipeline is provided under the MIT License.