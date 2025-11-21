rule create_multiqc_bam:
    input:
        "{species}/results/{reference}/analytics/",
    output:
        "{species}/results/{reference}/analytics/multiqc.html",
        directory("{species}/results/{reference}/analytics/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/{reference}/analytics/multiqc.log",
    wrapper:
        "v7.9.0/bio/multiqc"