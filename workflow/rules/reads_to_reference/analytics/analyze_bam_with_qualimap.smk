rule analyze_bam_with_qualimap:
    input:
        # BAM aligned, splicing-aware, to reference genome
        bam="{species}/processed/{reference}/mapped/{individual}_{reference}_final.bam"
    output:
        directory("{species}/results/{reference}/analytics/{individual}/qualimap"),
        "{species}/results/{reference}/analytics/{individual}/qualimap/qualimapReport.html"
    log:
        "{species}/results/{reference}/analytics/{individual}/{individual}_qualimap.log",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=4096,
    wrapper:
        "v7.6.0/bio/qualimap/bamqc"