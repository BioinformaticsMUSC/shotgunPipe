# Quality Control Rules

rule fastqc:
    """
    Run FastQC on raw reads
    """
    input:
        fastq=lambda wildcards: f"{config['data_dir']}/{wildcards.sample}_R{wildcards.read}.fastq.gz"
    output:
        html="results/fastqc/{sample}_R{read}_fastqc.html",
        zip="results/fastqc/{sample}_R{read}_fastqc.zip"
    params:
        outdir="results/fastqc"
    threads:
        config["fastqc"]["threads"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastqc {input.fastq} --outdir {params.outdir} --threads {threads}
        """

rule multiqc:
    """
    Aggregate QC reports with MultiQC
    """
    input:
        fastqc_reports=expand("results/fastqc/{sample}_R{read}_fastqc.zip", 
                             sample=SAMPLES, read=[1,2]),
        kneaddata_logs=expand("results/kneaddata/{sample}_kneaddata.log", 
                             sample=SAMPLES)
    output:
        report="results/multiqc_report.html"
    params:
        outdir="results",
        title=config["multiqc"]["title"],
        comment=config["multiqc"]["comment"]
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {input} --outdir {params.outdir} \
            --title "{params.title}" \
            --comment "{params.comment}" \
            --force
        """