# Annotation Rules - Functional Annotation of MAGs

rule prokka_annotation:
    """
    Annotate MAGs using Prokka
    """
    input:
        mag="results/binning/{sample}/bins/{sample}.bin.{bin_id}.fa"
    output:
        gff="results/annotation/{sample}/prokka/{bin_id}/{bin_id}.gff",
        faa="results/annotation/{sample}/prokka/{bin_id}/{bin_id}.faa",
        ffn="results/annotation/{sample}/prokka/{bin_id}/{bin_id}.ffn"
    params:
        outdir="results/annotation/{sample}/prokka/{bin_id}",
        prefix="{bin_id}",
        genus=config["prokka"]["genus"],
        species=config["prokka"]["species"]
    threads:
        config["prokka"]["threads"]
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        prokka {input.mag} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --genus {params.genus} \
            --species {params.species} \
            --cpus {threads} \
            --force
        """

rule dram_annotate:
    """
    Annotate MAGs using DRAM for comprehensive functional annotation
    """
    input:
        bins_dir="results/binning/{sample}/bins"
    output:
        annotations="results/annotation/{sample}/dram/annotations.tsv",
        summary="results/annotation/{sample}/dram/genome_summaries.tsv"
    params:
        outdir="results/annotation/{sample}/dram",
        min_contig_size=config["dram"]["min_contig_size"]
    threads:
        config["dram"]["threads"]
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        # Annotate genomes
        DRAM.py annotate \
            -i '{input.bins_dir}/*.fa' \
            -o {params.outdir} \
            --min_contig_size {params.min_contig_size} \
            --threads {threads}
        
        # Distill annotations
        DRAM.py distill \
            -i {params.outdir}/annotations.tsv \
            -o {params.outdir}/genome_summaries
        """

rule eggnog_mapper:
    """
    Functional annotation using eggNOG-mapper
    """
    input:
        proteins="results/annotation/{sample}/prokka/{bin_id}/{bin_id}.faa"
    output:
        annotations="results/annotation/{sample}/eggnog/{bin_id}.emapper.annotations"
    params:
        outdir="results/annotation/{sample}/eggnog",
        db_dir=config["eggnog"]["database_dir"],
        tax_scope=config["eggnog"]["tax_scope"]
    threads:
        config["eggnog"]["threads"]
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        emapper.py -i {input.proteins} \
            --output {wildcards.bin_id} \
            --output_dir {params.outdir} \
            --data_dir {params.db_dir} \
            --tax_scope {params.tax_scope} \
            --go_evidence non-electronic \
            --target_orthologs all \
            --seed_ortholog_evalue 0.001 \
            --seed_ortholog_score 60 \
            --override \
            --cpu {threads}
        """

rule kegg_pathway_analysis:
    """
    KEGG pathway analysis for annotated MAGs
    """
    input:
        annotations="results/annotation/{sample}/eggnog/{bin_id}.emapper.annotations"
    output:
        kegg_summary="results/annotation/{sample}/kegg/{bin_id}_kegg_summary.tsv"
    params:
        outdir="results/annotation/{sample}/kegg"
    conda:
        "../envs/annotation.yaml"
    script:
        "../scripts/kegg_analysis.py"

rule cazy_annotation:
    """
    Carbohydrate-Active enZyme (CAZy) annotation using dbCAN
    """
    input:
        proteins="results/annotation/{sample}/prokka/{bin_id}/{bin_id}.faa"
    output:
        cazy="results/annotation/{sample}/cazy/{bin_id}_cazy.out"
    params:
        db_dir=config["cazy"]["database_dir"]
    threads:
        config["cazy"]["threads"]
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        mkdir -p results/annotation/{wildcards.sample}/cazy
        
        run_dbcan {input.proteins} protein \
            --out_dir results/annotation/{wildcards.sample}/cazy/{wildcards.bin_id}_dbcan \
            --db_dir {params.db_dir} \
            --tools all \
            --dia_cpu {threads} \
            --hmm_cpu {threads}
        
        # Extract main results
        cp results/annotation/{wildcards.sample}/cazy/{wildcards.bin_id}_dbcan/overview.txt {output.cazy}
        """

rule combine_annotations:
    """
    Combine all annotation results for each sample
    """
    input:
        prokka_gff=lambda wildcards: expand("results/annotation/{sample}/prokka/{bin_id}/{bin_id}.gff", 
                                           sample=wildcards.sample, 
                                           bin_id=get_bin_ids(wildcards.sample)),
        dram_summary="results/annotation/{sample}/dram/genome_summaries.tsv"
    output:
        combined="results/annotation/{sample}/combined_annotations.tsv"
    script:
        "../scripts/combine_annotations.py"