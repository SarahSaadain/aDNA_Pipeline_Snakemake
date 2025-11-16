####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def filter_reads_by_quality_input_sample(wildcards):
    # Determine if the sample is paired-end or single-end
    reads = remove_adapters_type_with_fastp_input_sample(wildcards)
    if len(reads) == 2:
        # Paired-end: use the merged reads from fastp_pe
        return [f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.pe.fastq.gz"]
    else:
        # Single-end: use the trimmed reads from fastp_se
        return [f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.se.fastq.gz"]
 
####################################################
# Snakemake rules
####################################################

# Rule: Quality filtering of reads using fastp
rule filter_reads_by_quality:
    input:
        sample=filter_reads_by_quality_input_sample
    output:
        trimmed=temp("{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz"),
        failed=temp("{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.failed.fastq.gz"),
        html=report(
            "{species}/results/reads/reads_quality_filtered/fastp_report/{sample}_quality_filtered.html",
            #caption="../report/fastp.rst",
            category="quality control",
            subcategory="fastp",
        ),
        json="{species}/results/reads/reads_quality_filtered/fastp_report/{sample}_quality_filtered.json",
    message: "Quality filtering reads in {input.sample}"
    log:
        "{species}/logs/reads/reads_quality_filtered/{sample}_quality_filtered.log",
    params:
        extra=f"--disable_adapter_trimming --qualified_quality_phred {config['pipeline']['raw_reads_processing']['quality_filtering'].get('settings', {}).get('min_quality','15')} --length_required {config['pipeline']['raw_reads_processing']['quality_filtering'].get('settings', {}).get('min_length','15')} --unqualified_percent_limit 40 --n_base_limit 5"
    threads: 15
    wrapper:
        "v7.5.0/bio/fastp"