#
rule prepare_custom_data_breadth:
    input:
        csv="{species}/results/{reference}/analytics/{individual}/coverage/{individual}_{reference}_coverage_analysis.csv"
    output:
        tsv="{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_coverage_analysis.tsv"
    params:
        individual="{individual}"
    script:
        "../../scripts/processing_summary/prepare_custom_data_breadth.py"