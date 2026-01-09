# Statistical Analysis Rules

rule prepare_metadata:
    """
    Prepare metadata file for MaAsLin2 analysis
    """
    input:
        config["sample_sheet"]
    output:
        "results/maaslin2/metadata.tsv"
    shell:
        """
        # Convert CSV to TSV and ensure proper formatting
        sed 's/,/\t/g' {input} > {output}
        """

rule maaslin2_taxonomy:
    """
    Statistical analysis of taxonomic profiles using MaAsLin2
    """
    input:
        abundance="results/metaphlan/merged_abundance_table.txt",
        metadata="results/maaslin2/metadata.tsv"
    output:
        directory("results/maaslin2/taxonomy")
    params:
        min_abundance=config["maaslin2"]["min_abundance"],
        min_prevalence=config["maaslin2"]["min_prevalence"],
        normalization=config["maaslin2"]["normalization"],
        transform=config["maaslin2"]["transform"],
        analysis_method=config["maaslin2"]["analysis_method"],
        max_significance=config["maaslin2"]["max_significance"]
    conda:
        "../envs/maaslin2.yaml"
    script:
        "../scripts/run_maaslin2.R"

rule maaslin2_function:
    """
    Statistical analysis of functional profiles using MaAsLin2  
    """
    input:
        abundance="results/humann/merged_pathabundance_relab.tsv",
        metadata="results/maaslin2/metadata.tsv"
    output:
        directory("results/maaslin2/function")
    params:
        min_abundance=config["maaslin2"]["min_abundance"],
        min_prevalence=config["maaslin2"]["min_prevalence"], 
        normalization=config["maaslin2"]["normalization"],
        transform=config["maaslin2"]["transform"],
        analysis_method=config["maaslin2"]["analysis_method"],
        max_significance=config["maaslin2"]["max_significance"]
    conda:
        "../envs/maaslin2.yaml"
    script:
        "../scripts/run_maaslin2.R"