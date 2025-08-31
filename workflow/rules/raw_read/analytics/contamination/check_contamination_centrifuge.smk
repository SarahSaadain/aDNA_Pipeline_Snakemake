rule centrifuge_analysis:
    input:
        fastq = "{species}/processed/quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        centrifuge_out = "{species}/results/qualitycontrol/centrifuge/{sample}_centrifuge.txt",
        report = "{species}/results/qualitycontrol/centrifuge/{sample}_report.tsv"
    params:
        db = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"]["database"],
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"]["executable"]
    threads: workflow.cores
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"]["conda_env"]
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
        centrifuge_out = "{species}/processed/qualitycontrol/centrifuge/{sample}.centrifuge.txt"
    output:
        taxon_counts = "{species}/results/qualitycontrol/centrifuge/{sample}.taxon_counts.txt"
    shell:
        r"""
        awk '$3 != 0 {{print $3}}' {input.centrifuge_out} \
            | sort \
            | uniq -c \
            | sort -nr \
            > {output.taxon_counts}
        """