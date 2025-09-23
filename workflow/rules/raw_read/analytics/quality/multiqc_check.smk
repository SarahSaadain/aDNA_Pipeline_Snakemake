#import input_manager as sm

# Rule: Run MultiQC on raw FastQC outputs
rule multiqc_raw:
    input:
        lambda wc: get_input_multiqc_raw(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/raw/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/raw/multiqc.log",
    message: "Running MultiQC on raw FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

# Rule: Run MultiQC on trimmed FastQC outputs
rule multiqc_trimmed:
    input:
        lambda wc: get_input_multiqc_trimmed(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/trimmed/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/trimmed/multiqc.log",
    message: "Running MultiQC on trimmed FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

# Rule: Run MultiQC on quality-filtered FastQC outputs
rule multiqc_quality_filtered:
    input:
        lambda wc: get_input_multiqc_quality_filtered(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/quality_filtered/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/quality_filtered/multiqc.log",
    message: "Running MultiQC on quality-filtered FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

# Rule: Run MultiQC on merged FastQC outputs
rule multiqc_merged:
    input:
        lambda wc: get_input_multiqc_merged(wc.species)
    output:
        "{species}/results/qualitycontrol/multiqc/merged/multiqc.html",
        #directory("{species}/results/qualitycontrol/multiqc/raw//multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/qualitycontrol/multiqc/merged/multiqc.log",
    message: "Running MultiQC on merged FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"