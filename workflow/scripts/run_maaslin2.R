#!/usr/bin/env Rscript

# MaAsLin2 Statistical Analysis Script
# This script runs MaAsLin2 on abundance tables

library(Maaslin2)
library(optparse)

# Parse Snakemake inputs
input_data <- snakemake@input[["abundance"]]
input_metadata <- snakemake@input[["metadata"]]
output_dir <- snakemake@output[[1]]

# Parse parameters
min_abundance <- snakemake@params[["min_abundance"]]
min_prevalence <- snakemake@params[["min_prevalence"]]
normalization <- snakemake@params[["normalization"]]
transform <- snakemake@params[["transform"]]
analysis_method <- snakemake@params[["analysis_method"]]
max_significance <- snakemake@params[["max_significance"]]

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read data
cat("Reading abundance data from:", input_data, "\n")
abundance_data <- read.table(input_data, header = TRUE, sep = "\t", 
                           row.names = 1, check.names = FALSE)

cat("Reading metadata from:", input_metadata, "\n")
metadata <- read.table(input_metadata, header = TRUE, sep = "\t",
                      row.names = 1, check.names = FALSE)

# Transpose abundance data (MaAsLin2 expects samples as rows, features as columns)
abundance_data_t <- t(abundance_data)

# Match sample names between abundance and metadata
common_samples <- intersect(rownames(abundance_data_t), rownames(metadata))
abundance_data_matched <- abundance_data_t[common_samples, ]
metadata_matched <- metadata[common_samples, ]

cat("Number of samples after matching:", length(common_samples), "\n")
cat("Number of features:", ncol(abundance_data_matched), "\n")

# Run MaAsLin2
cat("Running MaAsLin2 analysis...\n")

# Determine fixed effects (all columns except sample identifier)
fixed_effects <- colnames(metadata_matched)

fit_data <- Maaslin2(
  input_data = abundance_data_matched,
  input_metadata = metadata_matched,
  output = output_dir,
  fixed_effects = fixed_effects,
  min_abundance = min_abundance,
  min_prevalence = min_prevalence,
  normalization = normalization,
  transform = transform,
  analysis_method = analysis_method,
  max_significance = max_significance,
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

cat("MaAsLin2 analysis completed. Results saved to:", output_dir, "\n")