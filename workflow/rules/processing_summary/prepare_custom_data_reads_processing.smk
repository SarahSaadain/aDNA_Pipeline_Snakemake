rule prepare_custom_data_reads_processing:
    input:
        reads = "{species}/results/reads/statistics/{species}_reads_counts.csv",
        dedup = "{species}/results/{reference}/analytics/{individual}/dedup/{individual}_{reference}_final.dedup.json"
    output:
        "{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_reads_processing_summary.tsv",
    params:
        individual="{individual}",
        reference="{reference}",
    script:
        "../../scripts/processing_summary/prepare_custom_data_reads_processing.py"
        