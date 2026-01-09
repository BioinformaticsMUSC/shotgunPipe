#!/usr/bin/env python3

"""
PCA Analysis Script for Taxonomic Data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

def main():
    # Read input files
    abundance_file = snakemake.input.abundance_table
    metadata_file = snakemake.input.metadata
    
    # Output files
    pca_plot_file = snakemake.output.pca_plot
    pca_data_file = snakemake.output.pca_data
    
    # Read abundance data
    abundance_df = pd.read_csv(abundance_file, sep='\t', index_col=0)
    
    # Read metadata
    metadata_df = pd.read_csv(metadata_file, index_col=0)
    
    # Filter to species level only (remove higher taxonomic levels)
    species_rows = abundance_df.index.str.contains('s__') & ~abundance_df.index.str.contains('t__')
    abundance_species = abundance_df[species_rows]
    
    # Transpose so samples are rows and species are columns
    abundance_t = abundance_species.T
    
    # Match samples between abundance and metadata
    common_samples = list(set(abundance_t.index) & set(metadata_df.index))
    abundance_matched = abundance_t.loc[common_samples]
    metadata_matched = metadata_df.loc[common_samples]
    
    # Remove species with zero abundance across all samples
    abundance_filtered = abundance_matched.loc[:, (abundance_matched != 0).any(axis=0)]
    
    # Apply log transformation (add pseudocount to avoid log(0))
    abundance_log = np.log10(abundance_filtered + 1e-6)
    
    # Standardize the data
    scaler = StandardScaler()
    abundance_scaled = scaler.fit_transform(abundance_log)
    
    # Perform PCA
    pca = PCA(n_components=min(abundance_scaled.shape[0]-1, abundance_scaled.shape[1], 10))
    pca_result = pca.fit_transform(abundance_scaled)
    
    # Create PCA dataframe
    pca_df = pd.DataFrame(pca_result, 
                         columns=[f'PC{i+1}' for i in range(pca_result.shape[1])],
                         index=abundance_matched.index)
    
    # Add metadata to PCA results
    pca_with_metadata = pca_df.join(metadata_matched)
    
    # Save PCA data
    pca_with_metadata.to_csv(pca_data_file, sep='\t')
    
    # Create PCA plot
    plt.figure(figsize=(12, 8))
    
    # Check if we have condition information for coloring
    if 'condition' in metadata_matched.columns:
        conditions = metadata_matched['condition'].unique()
        colors = plt.cm.Set1(np.linspace(0, 1, len(conditions)))
        
        for i, condition in enumerate(conditions):
            mask = metadata_matched['condition'] == condition
            plt.scatter(pca_result[mask, 0], pca_result[mask, 1], 
                       c=[colors[i]], label=condition, s=100, alpha=0.7)
        plt.legend()
    else:
        plt.scatter(pca_result[:, 0], pca_result[:, 1], s=100, alpha=0.7)
    
    # Add sample labels
    for i, sample in enumerate(abundance_matched.index):
        plt.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    plt.title('PCA of Taxonomic Profiles')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(pca_plot_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"PCA analysis completed. Explained variance by top 5 PCs:")
    for i in range(min(5, len(pca.explained_variance_ratio_))):
        print(f"  PC{i+1}: {pca.explained_variance_ratio_[i]:.1%}")

if __name__ == "__main__":
    main()