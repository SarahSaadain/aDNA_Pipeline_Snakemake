#
rule summarize_coverage:
    input:
        "{species}/results/{reference}/analytics/{individual}/coverage/{individual}_{reference}_coverage_analysis.csv"
    output:
        "{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_coverage_summary.tsv"
    params:
        individual="{individual}",  
        reference="{reference}"
    script:
        "../../scripts/processing_summary/summarize_coverage.py"

rule prepare_custom_data_depth:
    input:
        csv="{species}/results/{reference}/analytics/{individual}/coverage/{individual}_{reference}_coverage_analysis.csv"
    output:
        avg="{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_depth_coverage_avg.csv",
        max="{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_depth_coverage_max.csv"
    script:
        "../../scripts/processing_summary/prepare_custom_data_depth.py"