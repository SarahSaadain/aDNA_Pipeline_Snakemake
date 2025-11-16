####################################################
# Snakemake rules
####################################################

# Rule: Run Kraken2 for contamination analysis
rule analyze_contamination_with_kraken:
    input:
        fastq = "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        kraken_out = "{species}/processed/contamination_analysis/kraken/{sample}_kraken.tsv"
    params:
        db = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["kraken"]["settings"]["database"],
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["kraken"]["settings"]["executable"]
    threads: 15
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

