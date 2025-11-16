####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def run_fastqc_adapter_removed_input(wildcards):
    # Determine if the sample is paired-end or single-end
    reads = remove_adapters_type_with_fastp_input_sample(wildcards)
    if len(reads) == 2:
        # Paired-end: use the merged reads from fastp_pe
        return f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.pe.fastq.gz"
    else:
        # Single-end: use the trimmed reads from fastp_se
        return f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.se.fastq.gz"

####################################################
# Snakemake rules
####################################################

# Rule: Run FastQC on raw reads
rule run_fastqc_raw:
    input:
        "{species}/raw/reads/{sample}.fastq.gz"
    output:
        html="{species}/results/reads/reads_raw/fastqc/{sample}_fastqc.html",
        zip="{species}/results/reads/reads_raw/fastqc/{sample}_fastqc.zip"
    message: "Running FastQC on raw reads for sample {wildcards.sample} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/logs/reads/reads_raw/fastqc/{sample}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.5.0/bio/fastqc"

# Rule: Run FastQC on adapter-trimmed reads
rule run_fastqc_adapter_removed:
    input:
        run_fastqc_adapter_removed_input
    output:
        html="{species}/results/reads/reads_trimmed/fastqc/{sample}_trimmed_fastqc.html",
        zip="{species}/results/reads/reads_trimmed/fastqc/{sample}_trimmed_fastqc.zip"
    message: "Running FastQC on adapter-trimmed reads for sample {wildcards.sample} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/logs/reads/reads_trimmed/fastqc/{sample}_trimmed_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.5.0/bio/fastqc"

# Rule: Run FastQC on quality-filtered reads
rule run_fastqc_quality_filtered:
    input:
        "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz"
    output:
        html="{species}/results/reads/reads_quality_filtered/fastqc/{sample}_quality_filtered_fastqc.html",
        zip="{species}/results/reads/reads_quality_filtered/fastqc/{sample}_quality_filtered_fastqc.zip"
    message: "Running FastQC on quality-filtered reads for sample {wildcards.sample} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/logs/reads/reads_quality_filtered/fastqc/{sample}_quality_filtered_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.5.0/bio/fastqc"

# Rule: Run FastQC on merged reads
rule run_fastqc_merged:
    input:
        merged="{species}/processed/reads/reads_merged/{individual}.fastq.gz"
    output:
        html="{species}/results/reads/reads_merged/fastqc/{individual}_fastqc.html",
        zip="{species}/results/reads/reads_merged/fastqc/{individual}_fastqc.zip"
    message: "Running FastQC on merged reads for individual {wildcards.individual} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/logs/reads/reads_merged/fastqc/{individual}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.2.0/bio/fastqc"