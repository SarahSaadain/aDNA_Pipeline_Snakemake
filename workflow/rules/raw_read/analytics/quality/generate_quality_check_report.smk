#import scripts.setup_snakemake as setup
#import input_manager as im

# Rule: Generate quality check report for species
rule generate_qc_report:
    input:
        fastqc_raw = lambda wc: get_input_multiqc_raw(wc.species),
        fastqc_trimmed = lambda wc: get_input_multiqc_trimmed(wc.species),
        fastqc_quality_filtered = lambda wc: get_input_multiqc_quality_filtered(wc.species),
        multiqc_raw = "{species}/results/qualitycontrol/multiqc/raw/multiqc.html",
        multiqc_trimmed = "{species}/results/qualitycontrol/multiqc/raw/multiqc.html",
        multiqc_quality_filtered = "{species}/results/qualitycontrol/multiqc/merged/multiqc.html"
    output:
        report("{species}/results/qualitycontrol/quality_check_report_{species}.html")
    params:
        species="{species}"
    message: "Generating quality check report for species {params.species}"
    script:
        "../../../../scripts/raw_reads/analytics/quality/generate_quality_check_report.py"
