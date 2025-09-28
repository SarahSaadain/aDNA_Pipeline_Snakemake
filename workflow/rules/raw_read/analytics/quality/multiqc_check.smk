#import input_manager as sm

def get_input_multiqc_raw(wildcards):
    return get_expected_output_fastqc_raw(wildcards.species)

# Rule: Run MultiQC on raw FastQC outputs
rule multiqc_raw:
    input:
        get_input_multiqc_raw
    output:
        report("{species}/results/reads/{species}_multiqc_raw.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/logs/reads/{species}_multiqc_raw.log"
    message: "Running MultiQC on raw FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

def get_input_multiqc_trimmed(wildcards):
    return get_expected_output_fastqc_trimmed(wildcards.species)

# Rule: Run MultiQC on trimmed FastQC outputs
rule multiqc_trimmed:
    input:
        get_input_multiqc_trimmed
    output:
        report("{species}/results/reads/{species}_multiqc_trimmed.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/logs/reads/{species}_multiqc_trimmed.log"
    message: "Running MultiQC on trimmed FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

def get_input_multiqc_quality_filtered(wildcards):
    return get_expected_output_fastqc_quality_filtered(wildcards.species)

# Rule: Run MultiQC on quality-filtered FastQC outputs
rule multiqc_quality_filtered:
    input:
        get_input_multiqc_quality_filtered
    output:
        report("{species}/results/reads/{species}_multiqc_quality_filtered.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/logs/reads/{species}_multiqc_quality_filtered.log"
    message: "Running MultiQC on quality-filtered FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"

def get_input_multiqc_merged(wildcards):
    return get_expected_output_fastqc_merged(wildcards.species)

# Rule: Run MultiQC on merged FastQC outputs
rule multiqc_merged:
    input:
        get_input_multiqc_merged
    output:
        report("{species}/results/reads/{species}_multiqc_merged.html",
            category="quality control",
            subcategory="multiqc",
        ),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/logs/reads/{species}_multiqc_merged.log"
    message: "Running MultiQC on merged FastQC outputs for species {wildcards.species}"
    wrapper:
        "v7.2.0/bio/multiqc"