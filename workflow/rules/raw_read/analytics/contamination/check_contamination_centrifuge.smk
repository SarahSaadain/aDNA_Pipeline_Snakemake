####################################################
# Snakemake rules
####################################################

rule analyze_contamination_with_centrifuge:
    input:
        fastq = "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        centrifuge_out = "{species}/results/contamination_analysis/centrifuge/{sample}_centrifuge.txt",
        report = "{species}/results/contamination_analysis/centrifuge/{sample}_report.tsv"
    params:
        db = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"]["database"],
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"]["executable"]
    threads: 15
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"]["conda_env"]
    message: "Running Centrifuge contamination analysis for {input.fastq}"
    shell:
        """
        # Run centrifuge
        {params.executable] \
            -x {params.db} \
            -U {input.fastq} \
            -S {output.centrifuge_out} \
            --report-file {output.report} \
            --threads {threads} \
            --seed 999
        """

rule taxon_counts:
    input:
        centrifuge_out = "{species}/results/contamination_analysis/centrifuge/{sample}.centrifuge.txt"
    output:
        taxon_counts = "{species}/results/contamination_analysis/centrifuge/{sample}.taxon_counts.txt"
    message: "Counting taxon occurrences in {input.centrifuge_out}"
    shell:
        r"""
        awk '$3 != 0 {{print $3}}' {input.centrifuge_out} \
            | sort \
            | uniq -c \
            | sort -nr \
            > {output.taxon_counts}
        """