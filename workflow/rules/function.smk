# Functional Profiling Rules

rule humann:
    """
    Functional profiling with HUMAnN3
    """
    input:
        reads="results/kneaddata/{sample}_concat.fastq",
        metaphlan_profile="results/metaphlan/{sample}_metaphlan.txt"
    output:
        genefamilies="results/humann/{sample}_genefamilies.tsv",
        pathabundance="results/humann/{sample}_pathabundance.tsv",
        pathcoverage="results/humann/{sample}_pathcoverage.tsv"
    params:
        outdir="results/humann",
        chocophlan_db=config["databases"]["chocophlan_db"],
        uniref_db=config["databases"]["uniref_db"],
        search_mode=config["humann"]["search_mode"],
        memory_use=config["humann"]["memory_use"],
        bypass_translated=config["humann"]["bypass_translated_search"]
    threads:
        config["humann"]["threads"]
    conda:
        "../envs/humann.yaml"
    shell:
        """
        humann --input {input.reads} \
            --output {params.outdir} \
            --nucleotide-database {params.chocophlan_db} \
            --protein-database {params.uniref_db} \
            --taxonomic-profile {input.metaphlan_profile} \
            --search-mode {params.search_mode} \
            --memory-use {params.memory_use} \
            --threads {threads} \
            {"--bypass-translated-search" if params.bypass_translated else ""} \
            --remove-temp-output
        """

rule humann_rename_tables:
    """
    Rename HUMAnN output files to include sample name
    """
    input:
        genefamilies="results/humann/{sample}_genefamilies.tsv",
        pathabundance="results/humann/{sample}_pathabundance.tsv", 
        pathcoverage="results/humann/{sample}_pathcoverage.tsv"
    output:
        genefamilies_renamed="results/humann/{sample}_genefamilies_renamed.tsv",
        pathabundance_renamed="results/humann/{sample}_pathabundance_renamed.tsv",
        pathcoverage_renamed="results/humann/{sample}_pathcoverage_renamed.tsv"
    conda:
        "../envs/humann.yaml"
    shell:
        """
        humann_rename_table --input {input.genefamilies} --output {output.genefamilies_renamed} --names uniref90
        humann_rename_table --input {input.pathabundance} --output {output.pathabundance_renamed} --names metacyc
        humann_rename_table --input {input.pathcoverage} --output {output.pathcoverage_renamed} --names metacyc
        """

rule merge_humann_genefamilies:
    """
    Merge HUMAnN gene families tables
    """
    input:
        genefamilies=expand("results/humann/{sample}_genefamilies.tsv", sample=SAMPLES)
    output:
        merged="results/humann/merged_genefamilies.tsv"
    conda:
        "../envs/humann.yaml"
    shell:
        """
        humann_join_tables --input results/humann --output {output.merged} --file_name genefamilies
        """

rule merge_humann_pathabundance:
    """
    Merge HUMAnN pathway abundance tables
    """
    input:
        pathabundance=expand("results/humann/{sample}_pathabundance.tsv", sample=SAMPLES)
    output:
        merged="results/humann/merged_pathabundance.tsv"
    conda:
        "../envs/humann.yaml"
    shell:
        """
        humann_join_tables --input results/humann --output {output.merged} --file_name pathabundance
        """

rule normalize_humann_tables:
    """
    Normalize HUMAnN tables to relative abundance
    """
    input:
        genefamilies="results/humann/merged_genefamilies.tsv",
        pathabundance="results/humann/merged_pathabundance.tsv"
    output:
        genefamilies_norm="results/humann/merged_genefamilies_relab.tsv",
        pathabundance_norm="results/humann/merged_pathabundance_relab.tsv"
    conda:
        "../envs/humann.yaml"
    shell:
        """
        humann_renorm_table --input {input.genefamilies} --output {output.genefamilies_norm} --units relab
        humann_renorm_table --input {input.pathabundance} --output {output.pathabundance_norm} --units relab
        """