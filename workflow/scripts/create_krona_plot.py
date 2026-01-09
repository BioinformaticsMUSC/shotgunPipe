#!/usr/bin/env python3

"""
Create Krona plot for taxonomic visualization
"""

import pandas as pd
import subprocess
import tempfile
import os

def main():
    # Input and output files
    metaphlan_file = snakemake.input.metaphlan
    krona_output = snakemake.output.krona
    
    # Read MetaPhlAn output
    df = pd.read_csv(metaphlan_file, sep='\t', comment='#', header=0)
    
    # Get the abundance column (usually the second column)
    abundance_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]
    
    # Create temporary file for Krona input
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_file:
        tmp_filename = tmp_file.name
        
        for idx, row in df.iterrows():
            taxonomy = row.iloc[0]  # First column is taxonomy
            abundance = row[abundance_col]  # Abundance column
            
            # Skip unclassified entries
            if 'unclassified' in taxonomy.lower():
                continue
                
            # Parse taxonomy string
            tax_levels = taxonomy.split('|')
            
            # Remove prefixes (k__, p__, c__, o__, f__, g__, s__)
            clean_levels = []
            for level in tax_levels:
                if '__' in level:
                    clean_levels.append(level.split('__')[1])
                else:
                    clean_levels.append(level)
            
            # Write to Krona input format: abundance followed by taxonomy levels
            if clean_levels and abundance > 0:
                tmp_file.write(f"{abundance}\t" + '\t'.join(clean_levels) + '\n')
    
    # Create Krona plot
    cmd = [
        'ktImportText',
        tmp_filename,
        '-o', krona_output
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Krona plot created successfully: {krona_output}")
    except subprocess.CalledProcessError as e:
        print(f"Error creating Krona plot: {e.stderr}")
        # Create a simple HTML file if Krona fails
        with open(krona_output, 'w') as f:
            f.write("<html><body><h1>Krona plot generation failed</h1></body></html>")
    finally:
        # Clean up temporary file
        os.unlink(tmp_filename)

if __name__ == "__main__":
    main()