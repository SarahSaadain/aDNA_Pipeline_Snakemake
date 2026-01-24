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

rule analyze_centrifuge_report_taxon_counts:
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

rule analyze_centrifuge_report_proportions:
    input:
        report = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_centrifuge_report.tsv",
        count_reads = "{species}/processed/reads/statistics/{sample}_quality_filtered.count"
    output:
        "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_proportions.tsv",
    params:
        sample = "{sample}"
    script:
        "../../../../scripts/raw_reads/analytics/contamination/check_contamination_ecmsd_script_analyze_centrifuge_report_proportions.py"

rule analyze_centrifuge_report_top_taxa:
    input:
        report = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_centrifuge_report.tsv",
    output:
        top10_unique = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_top10_unique_taxa.tsv",
        top10_total = "{species}/results/contamination_analysis/centrifuge/{individual}/{sample}/{sample}_top10_total_taxa.tsv"
    params:
        sample = "{sample}"
    script:
        "../../../../scripts/raw_reads/analytics/contamination/check_contamination_ecmsd_script_analyze_centrifuge_report_top_taxa.py"