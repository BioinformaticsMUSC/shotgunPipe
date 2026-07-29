# Database setup rules

rule setup_databases:
    """
    Download and prepare all reference databases used by the pipeline.
    """
    input:
        kneaddata_host="resources/databases/kneaddata/human/.complete",
        metaphlan_db="resources/databases/metaphlan/.complete",
        humann_chocophlan="resources/databases/humann/chocophlan/.complete",
        humann_uniref="resources/databases/humann/uniref90/.complete",
        eggnog_db="resources/databases/eggnog/.complete",
        dbcan_db="resources/databases/dbcan/.complete"


rule download_kneaddata_host_db:
    """
    Download the host genome database for KneadData.
    """
    output:
        complete="resources/databases/kneaddata/human/.complete"
    params:
        db_dir=config["databases"]["host_genome"]
    conda:
        "../envs/kneaddata.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        kneaddata_database --download human_genome bowtie2 {params.db_dir}
        touch {output.complete}
        """


rule download_metaphlan_db:
    """
    Download the MetaPhlAn database.
    """
    output:
        complete="resources/databases/metaphlan/.complete"
    params:
        db_dir=config["databases"]["metaphlan_db"]
    conda:
        "../envs/metaphlan.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        metaphlan --install --bowtie2db {params.db_dir}
        touch {output.complete}
        """


rule download_humann_chocophlan_db:
    """
    Download the HUMAnN ChocoPhlAn nucleotide database.
    """
    output:
        complete="resources/databases/humann/chocophlan/.complete"
    params:
        db_dir=config["databases"]["chocophlan_db"]
    conda:
        "../envs/humann.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        humann_databases --download chocophlan full {params.db_dir}
        touch {output.complete}
        """


rule download_humann_uniref_db:
    """
    Download the HUMAnN UniRef90 translated search database.
    """
    output:
        complete="resources/databases/humann/uniref90/.complete"
    params:
        db_dir=config["databases"]["uniref_db"]
    conda:
        "../envs/humann.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        humann_databases --download uniref uniref90_diamond {params.db_dir}
        touch {output.complete}
        """


rule download_eggnog_db:
    """
    Download the eggNOG database for eggNOG-mapper.
    """
    output:
        complete="resources/databases/eggnog/.complete"
    params:
        db_dir=config["eggnog"]["database_dir"]
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        download_eggnog_data.py -y --data_dir {params.db_dir}
        touch {output.complete}
        """


rule download_dbcan_db:
    """
    Download the dbCAN database used by run_dbcan.
    """
    output:
        complete="resources/databases/dbcan/.complete"
    params:
        db_dir=config["cazy"]["database_dir"]
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        mkdir -p {params.db_dir}
        run_dbcan database --db_dir {params.db_dir} --aws_s3
        touch {output.complete}
        """