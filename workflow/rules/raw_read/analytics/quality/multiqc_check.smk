import input_manager as sm

rule multiqc_raw:
    input:
        lambda wc: sm.get_input_multiqc_raw(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/raw/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/raw/multiqc.log",
    wrapper:
        "v7.2.0/bio/multiqc"

rule multiqc_trimmed:
    input:
        lambda wc: sm.get_input_multiqc_trimmed(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/trimmed/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/trimmed/multiqc.log",
    wrapper:
        "v7.2.0/bio/multiqc"

rule multiqc_quality_filtered:
    input:
        lambda wc: sm.get_input_multiqc_quality_filtered(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/quality_filtered/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/quality_filtered/multiqc.log",
    wrapper:
        "v7.2.0/bio/multiqc"

rule multiqc_merged:
    input:
        lambda wc: sm.get_input_multiqc_merged(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/merged/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/merged/multiqc.log",
    wrapper:
        "v7.2.0/bio/multiqc"