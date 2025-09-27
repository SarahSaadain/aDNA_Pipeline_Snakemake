#import input_manager as sm

# Rule: Run MultiQC on raw FastQC outputs
rule multiqc_raw:
    input:
        lambda wc: get_input_multiqc_raw(wc.species)
    output:
        report("{species}/results/reads/{species}_multiqc_raw.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/reads/{species}_multiqc_raw.log"
    message: "Running MultiQC on raw FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

# Rule: Run MultiQC on trimmed FastQC outputs
rule multiqc_trimmed:
    input:
        lambda wc: get_input_multiqc_trimmed(wc.species)
    output:
        report("{species}/results/reads/{species}_multiqc_trimmed.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/reads/{species}_multiqc_trimmed.log"
    message: "Running MultiQC on trimmed FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

# Rule: Run MultiQC on quality-filtered FastQC outputs
rule multiqc_quality_filtered:
    input:
        lambda wc: get_input_multiqc_quality_filtered(wc.species)
    output:
        report("{species}/results/reads/{species}_multiqc_quality_filtered.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/reads/{species}_multiqc_quality_filtered.log"
    message: "Running MultiQC on quality-filtered FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

# Rule: Run MultiQC on merged FastQC outputs
rule multiqc_merged:
    input:
        lambda wc: get_input_multiqc_merged(wc.species)
    output:
        report("{species}/results/reads/{species}_multiqc_merged.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/reads/{species}_multiqc_merged.log"
    message: "Running MultiQC on merged FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"