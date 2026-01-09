# Visualization Rules - Data Visualization and Exploration

rule krona_taxonomy:
    """
    Create Krona plot for taxonomic visualization
    """
    input:
        metaphlan="results/metaphlan/{sample}_metaphlan.txt"
    output:
        krona="results/visualization/{sample}_taxonomy_krona.html"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/create_krona_plot.py"

rule krona_function:
    """
    Create Krona plot for functional pathways
    """
    input:
        humann="results/humann/{sample}_pathabundance.tsv"
    output:
        krona="results/visualization/{sample}_function_krona.html"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/create_krona_function.py"

rule pca_analysis:
    """
    Principal Component Analysis of taxonomic data
    """
    input:
        abundance_table="results/metaphlan/merged_abundance_table.txt",
        metadata=config["sample_sheet"]
    output:
        pca_plot="results/visualization/pca_taxonomy.png",
        pca_data="results/visualization/pca_taxonomy_data.tsv"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/pca_analysis.py"

rule ordination_analysis:
    """
    NMDS and other ordination methods
    """
    input:
        abundance_table="results/metaphlan/merged_abundance_table.txt",
        metadata=config["sample_sheet"]
    output:
        nmds_plot="results/visualization/nmds_taxonomy.png",
        pcoa_plot="results/visualization/pcoa_taxonomy.png"
    params:
        distance_method=config["visualization"]["distance_method"]
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/ordination_analysis.R"

rule diversity_analysis:
    """
    Alpha and beta diversity analysis
    """
    input:
        abundance_table="results/metaphlan/merged_abundance_table.txt",
        metadata=config["sample_sheet"]
    output:
        alpha_diversity="results/visualization/alpha_diversity.png",
        beta_diversity="results/visualization/beta_diversity.png",
        diversity_stats="results/visualization/diversity_statistics.tsv"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/diversity_analysis.R"

rule taxonomy_barplot:
    """
    Create stacked barplot of taxonomic composition
    """
    input:
        abundance_table="results/metaphlan/merged_abundance_table.txt",
        metadata=config["sample_sheet"]
    output:
        barplot="results/visualization/taxonomy_barplot.png",
        top_taxa="results/visualization/top_taxa_summary.tsv"
    params:
        top_n=config["visualization"]["top_taxa"]
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/taxonomy_barplot.py"

rule functional_barplot:
    """
    Create stacked barplot of functional pathways
    """
    input:
        pathway_table="results/humann/merged_pathabundance_relab.tsv",
        metadata=config["sample_sheet"]
    output:
        barplot="results/visualization/pathway_barplot.png",
        top_pathways="results/visualization/top_pathways_summary.tsv"
    params:
        top_n=config["visualization"]["top_pathways"]
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/functional_barplot.py"

rule heatmap_analysis:
    """
    Create heatmaps for differential features
    """
    input:
        maaslin_taxonomy="results/maaslin2/taxonomy/significant_results.tsv",
        maaslin_function="results/maaslin2/function/significant_results.tsv",
        taxonomy_table="results/metaphlan/merged_abundance_table.txt",
        function_table="results/humann/merged_pathabundance_relab.tsv"
    output:
        taxonomy_heatmap="results/visualization/significant_taxa_heatmap.png",
        function_heatmap="results/visualization/significant_pathways_heatmap.png"
    params:
        max_features=config["visualization"]["max_heatmap_features"],
        cluster_samples=config["visualization"]["cluster_samples"]
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/create_heatmaps.py"

rule mag_quality_plot:
    """
    Visualize MAG quality metrics
    """
    input:
        checkm_results=expand("results/binning/{sample}/checkm/quality_summary.tsv", sample=SAMPLES)
    output:
        quality_plot="results/visualization/mag_quality.png",
        quality_summary="results/visualization/mag_quality_summary.tsv"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/mag_quality_viz.py"

rule assembly_stats_plot:
    """
    Visualize assembly statistics
    """
    input:
        assembly_stats=expand("results/assembly/{sample}/quast_stats", sample=SAMPLES)
    output:
        stats_plot="results/visualization/assembly_stats.png"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/assembly_stats_viz.py"

rule create_dashboard:
    """
    Create interactive HTML dashboard with all visualizations
    """
    input:
        pca_plot="results/visualization/pca_taxonomy.png",
        taxonomy_barplot="results/visualization/taxonomy_barplot.png",
        functional_barplot="results/visualization/pathway_barplot.png",
        alpha_diversity="results/visualization/alpha_diversity.png",
        taxonomy_heatmap="results/visualization/significant_taxa_heatmap.png",
        mag_quality="results/visualization/mag_quality.png",
        multiqc_report="results/multiqc_report.html"
    output:
        dashboard="results/visualization/metagenomics_dashboard.html"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/create_dashboard.py"