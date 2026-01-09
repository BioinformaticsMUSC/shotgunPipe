# Binning Rules - Metagenome-Assembled Genomes (MAGs)

rule metabat2_depth:
    """
    Calculate contig depths for MetaBAT2 binning
    """
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        r1="results/kneaddata/{sample}_kneaddata_paired_1.fastq",
        r2="results/kneaddata/{sample}_kneaddata_paired_2.fastq"
    output:
        bam="results/binning/{sample}/{sample}.sorted.bam",
        depth="results/binning/{sample}/{sample}.depth.txt"
    threads:
        config["metabat2"]["threads"]
    conda:
        "../envs/binning.yaml"
    shell:
        """
        # Index contigs
        bowtie2-build {input.contigs} results/binning/{wildcards.sample}/{wildcards.sample}_contigs
        
        # Map reads to contigs
        bowtie2 -x results/binning/{wildcards.sample}/{wildcards.sample}_contigs \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} | \
        samtools sort -@ {threads} -o {output.bam}
        
        # Index BAM file
        samtools index {output.bam}
        
        # Calculate depth
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam}
        """

rule metabat2_binning:
    """
    Bin contigs into MAGs using MetaBAT2
    """
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        depth="results/binning/{sample}/{sample}.depth.txt"
    output:
        bins_dir=directory("results/binning/{sample}/bins")
    params:
        outprefix="results/binning/{sample}/bins/{sample}.bin",
        min_contig=config["metabat2"]["min_contig_len"],
        min_cv=config["metabat2"]["min_cv"],
        min_cv_sum=config["metabat2"]["min_cv_sum"]
    threads:
        config["metabat2"]["threads"]
    conda:
        "../envs/binning.yaml"
    shell:
        """
        mkdir -p {output.bins_dir}
        
        metabat2 -i {input.contigs} -a {input.depth} \
            -o {params.outprefix} \
            --minContig {params.min_contig} \
            --minCV {params.min_cv} \
            --minCVSum {params.min_cv_sum} \
            --numThreads {threads}
        """

rule concoct_binning:
    """
    Alternative binning using CONCOCT
    """
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        bam="results/binning/{sample}/{sample}.sorted.bam"
    output:
        bins_dir=directory("results/binning/{sample}/concoct_bins"),
        coverage="results/binning/{sample}/concoct_coverage.tsv"
    params:
        chunk_size=config["concoct"]["chunk_size"],
        overlap_size=config["concoct"]["overlap_size"]
    threads:
        config["concoct"]["threads"]
    conda:
        "../envs/binning.yaml"
    shell:
        """
        mkdir -p {output.bins_dir}
        
        # Cut contigs into smaller pieces
        cut_up_fasta.py {input.contigs} \
            -c {params.chunk_size} \
            -o {params.overlap_size} \
            --merge_last \
            -b results/binning/{wildcards.sample}/contigs_10K.bed \
            > results/binning/{wildcards.sample}/contigs_10K.fa
        
        # Generate coverage table
        concoct_coverage_table.py results/binning/{wildcards.sample}/contigs_10K.bed \
            {input.bam} > {output.coverage}
        
        # Run CONCOCT
        concoct --composition_file results/binning/{wildcards.sample}/contigs_10K.fa \
            --coverage_file {output.coverage} \
            -b {output.bins_dir}/ \
            --threads {threads}
        
        # Merge subcontig clustering into original contig clustering
        merge_cutup_clustering.py results/binning/{wildcards.sample}/clustering_gt1000.csv \
            > results/binning/{wildcards.sample}/clustering_merged.csv
        
        # Extract bins as fasta files
        extract_fasta_bins.py {input.contigs} \
            results/binning/{wildcards.sample}/clustering_merged.csv \
            --output_path {output.bins_dir}
        """

rule checkm_quality:
    """
    Assess MAG quality using CheckM
    """
    input:
        bins_dir="results/binning/{sample}/bins"
    output:
        quality="results/binning/{sample}/checkm/quality_summary.tsv"
    params:
        checkm_dir="results/binning/{sample}/checkm",
        extension="fa"
    threads:
        config["checkm"]["threads"]
    conda:
        "../envs/binning.yaml"
    shell:
        """
        mkdir -p {params.checkm_dir}
        
        checkm lineage_wf -f {output.quality} \
            -t {threads} \
            -x {params.extension} \
            {input.bins_dir} \
            {params.checkm_dir}
        """

rule dastool_binning:
    """
    Improve binning results using DAS Tool
    """
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        metabat_bins="results/binning/{sample}/bins",
        concoct_bins="results/binning/{sample}/concoct_bins"
    output:
        scaffolds2bin="results/binning/{sample}/dastool_scaffolds2bin.tsv",
        optimized_bins=directory("results/binning/{sample}/dastool_bins")
    params:
        outprefix="results/binning/{sample}/dastool"
    threads:
        config["dastool"]["threads"]
    conda:
        "../envs/binning.yaml"
    shell:
        """
        # Prepare scaffolds2bin files
        Fasta_to_Scaffolds2Bin.sh -i {input.metabat_bins} \
            -e fa > {params.outprefix}_metabat.scaffolds2bin.tsv
        Fasta_to_Scaffolds2Bin.sh -i {input.concoct_bins} \
            -e fa > {params.outprefix}_concoct.scaffolds2bin.tsv
        
        # Run DAS Tool
        DAS_Tool -i {params.outprefix}_metabat.scaffolds2bin.tsv,{params.outprefix}_concoct.scaffolds2bin.tsv \
            -l metabat,concoct \
            -c {input.contigs} \
            -o {params.outprefix} \
            --write_bins \
            --threads {threads}
        
        # Create output directory
        mkdir -p {output.optimized_bins}
        cp {params.outprefix}_DASTool_bins/*.fa {output.optimized_bins}/ || true
        """