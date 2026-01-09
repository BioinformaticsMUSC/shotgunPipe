# Assembly Rules - Metagenome Assembly

rule megahit_assembly:
    """
    Assemble metagenomes using MEGAHIT
    """
    input:
        r1="results/kneaddata/{sample}_kneaddata_paired_1.fastq",
        r2="results/kneaddata/{sample}_kneaddata_paired_2.fastq"
    output:
        contigs="results/assembly/{sample}/final.contigs.fa",
        log="results/assembly/{sample}/log"
    params:
        outdir="results/assembly/{sample}",
        min_contig_len=config["megahit"]["min_contig_len"],
        k_list=config["megahit"]["k_list"],
        memory=config["megahit"]["memory"]
    threads:
        config["megahit"]["threads"]
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        rm -rf {params.outdir}
        megahit -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            --min-contig-len {params.min_contig_len} \
            --k-list {params.k_list} \
            --memory {params.memory} \
            --num-cpu-threads {threads} \
            --force
        """

rule metaspades_assembly:
    """
    Assemble metagenomes using metaSPAdes (alternative to MEGAHIT)
    """
    input:
        r1="results/kneaddata/{sample}_kneaddata_paired_1.fastq",
        r2="results/kneaddata/{sample}_kneaddata_paired_2.fastq"
    output:
        contigs="results/assembly_spades/{sample}/contigs.fasta",
        scaffolds="results/assembly_spades/{sample}/scaffolds.fasta"
    params:
        outdir="results/assembly_spades/{sample}",
        memory=config["metaspades"]["memory"],
        k_list=config["metaspades"]["k_list"]
    threads:
        config["metaspades"]["threads"]
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        rm -rf {params.outdir}
        spades.py --meta \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            -m {params.memory} \
            -t {threads} \
            -k {params.k_list}
        """

rule assembly_stats:
    """
    Calculate assembly statistics using QUAST
    """
    input:
        contigs="results/assembly/{sample}/final.contigs.fa"
    output:
        stats_dir=directory("results/assembly/{sample}/quast_stats")
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        quast.py {input.contigs} -o {output.stats_dir} --threads {threads}
        """

rule combine_assemblies:
    """
    Combine all individual assemblies for downstream analysis
    """
    input:
        assemblies=expand("results/assembly/{sample}/final.contigs.fa", sample=SAMPLES)
    output:
        combined="results/assembly/all_contigs_combined.fa"
    shell:
        """
        cat {input.assemblies} > {output.combined}
        """