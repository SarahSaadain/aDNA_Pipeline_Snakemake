# Rule: Run Kraken2 for contamination analysis
rule kraken_analysis:
    input:
        fastq = "{species}/processed/quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        kraken_out = "{species}/processed/qualitycontrol/kraken/{sample}_kraken.tsv"
    params:
        db = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["kraken"]["settings"]["database"],
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["kraken"]["settings"]["executable"]
    threads: workflow.cores
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["kraken"]["settings"]["conda_env"]
    message: "Running Kraken2 contamination analysis for {input.fastq}"
    shell:
        """
        # Run Kraken2
        {params.executable} \
            --db {params.db} \
            --threads {threads} \
            --gzip-compressed \
            --output {output.kraken_out} \
            {input.fastq}
        """

