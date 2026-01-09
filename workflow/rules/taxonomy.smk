# Taxonomic Profiling Rules

rule metaphlan:
    """
    Taxonomic profiling with MetaPhlAn4
    """
    input:
        reads="results/kneaddata/{sample}_concat.fastq"
    output:
        profile="results/metaphlan/{sample}_metaphlan.txt",
        bowtie2_out="results/metaphlan/{sample}_metaphlan.bowtie2.bz2"
    params:
        db=config["databases"]["metaphlan_db"],
        analysis_type=config["metaphlan"]["analysis_type"],
        add_viruses="--add_viruses" if config["metaphlan"]["add_viruses"] else "",
        unknown_est="--unknown_estimation" if config["metaphlan"]["unknown_estimation"] else ""
    threads:
        config["metaphlan"]["threads"]
    conda:
        "../envs/metaphlan.yaml"
    shell:
        """
        metaphlan {input.reads} \
            --input_type fastq \
            --bowtie2db {params.db} \
            --bowtie2out {output.bowtie2_out} \
            --nproc {threads} \
            --analysis_type {params.analysis_type} \
            {params.add_viruses} \
            {params.unknown_est} \
            --output_file {output.profile}
        """

rule merge_metaphlan:
    """
    Merge MetaPhlAn profiles into a single abundance table
    """
    input:
        profiles=expand("results/metaphlan/{sample}_metaphlan.txt", sample=SAMPLES)
    output:
        merged="results/metaphlan/merged_abundance_table.txt"
    conda:
        "../envs/metaphlan.yaml"
    shell:
        """
        merge_metaphlan_tables.py {input.profiles} > {output.merged}
        """