rule analyze_bam_with_picard_duplicates:
    input:
        bams="{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam",
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam=temp("{species}/results/{reference}/analytics/{individual}/picard_duplicates/{individual}_{reference}_picard_dedup.bam"),
        metrics="{species}/results/{reference}/analytics/{individual}/picard_duplicates/{individual}_{reference}_metrics.txt",
    log:
        "{species}/results/{reference}/analytics/{individual}/picard_duplicates/{individual}_{reference}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v7.6.0/bio/picard/markduplicates"