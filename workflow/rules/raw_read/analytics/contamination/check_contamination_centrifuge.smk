####################################################
# Snakemake rules
####################################################

rule analyze_contamination_with_centrifuge:
    input:
        fastq = "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        output = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_centrifuge_output.tsv",
        report = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_centrifuge_report.tsv"
    params:
        index = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"]["index"],
    threads: 15
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["centrifuge"]["settings"].get("conda_env", "../../../../envs/centrifuge.yaml")
    message: "Running Centrifuge contamination analysis for {input.fastq}"
    shell:
        """
        # Run centrifuge
        centrifuge \
            -x {params.index} \
            -U {input.fastq} \
            -S {output.output} \
            --report-file {output.report} \
            --threads {threads}
        """

rule taxon_counts:
    input:
        centrifuge_out = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_centrifuge_output.tsv"
    output:
        taxon_counts = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_taxon_counts.tsv"
    message: "Counting taxon occurrences in {input.centrifuge_out}"
    shell:
        r"""
        awk '$3 != 0 {{print $3}}' {input.centrifuge_out} \
            | sort \
            | uniq -c \
            | sort -nr \
            > {output.taxon_counts}
        """