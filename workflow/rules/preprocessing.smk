# Preprocessing Rules - Host Removal

rule kneaddata:
    """
    Remove host contamination and low-quality reads using KneadData
    """
    input:
        r1=lambda wildcards: get_fastq_path(wildcards.sample, 1),
        r2=lambda wildcards: get_fastq_path(wildcards.sample, 2)
    output:
        r1_clean="results/kneaddata/{sample}_kneaddata_paired_1.fastq",
        r2_clean="results/kneaddata/{sample}_kneaddata_paired_2.fastq",
        log="results/kneaddata/{sample}_kneaddata.log"
    params:
        outdir="results/kneaddata",
        host_db=config["databases"]["host_genome"],
        trimmomatic_opts=config["kneaddata"]["trimmomatic_options"],
        bowtie2_opts=config["kneaddata"]["bowtie2_options"]
    threads:
        config["kneaddata"]["threads"]
    conda:
        "../envs/kneaddata.yaml"
    shell:
        """
        kneaddata --input {input.r1} --input {input.r2} \
            --output {params.outdir} \
            --reference-db {params.host_db} \
            --threads {threads} \
            --processes {threads} \
            --trimmomatic-options "{params.trimmomatic_opts}" \
            --bowtie2-options "{params.bowtie2_opts}" \
            --log {output.log} \
            --remove-intermediate-output
        """

rule concatenate_reads:
    """
    Concatenate paired-end reads for downstream analysis
    """
    input:
        r1="results/kneaddata/{sample}_kneaddata_paired_1.fastq",
        r2="results/kneaddata/{sample}_kneaddata_paired_2.fastq"
    output:
        concat="results/kneaddata/{sample}_concat.fastq"
    shell:
        """
        cat {input.r1} {input.r2} > {output.concat}
        """