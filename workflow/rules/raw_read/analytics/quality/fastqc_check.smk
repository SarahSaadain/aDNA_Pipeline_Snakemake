def get_fastqc_adapter_removed_input_reads(wildcards):
    # Determine if the sample is paired-end or single-end
    reads = get_adapter_removal_input_reads(wildcards)
    if len(reads) == 2:
        # Paired-end: use the merged reads from fastp_pe
        return f"{wildcards.species}/processed/trimmed/{wildcards.sample}_trimmed.pe.fastq.gz"
    else:
        # Single-end: use the trimmed reads from fastp_se
        return f"{wildcards.species}/processed/trimmed/{wildcards.sample}_trimmed.se.fastq.gz"

# Rule: Run FastQC on raw reads
rule fastqc_raw:
    input:
        "{species}/raw/reads/{sample}.fastq.gz"
    output:
        html="{species}/results/qualitycontrol/fastqc/raw/{sample}_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/raw/{sample}_fastqc.zip"
    message: "Running FastQC on raw reads for sample {wildcards.sample} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/results/qualitycontrol/fastqc/raw/{sample}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.5.0/bio/fastqc"

# Rule: Run FastQC on adapter-trimmed reads
rule fastqc_adapter_removed:
    input:
        get_fastqc_adapter_removed_input_reads
    output:
        html="{species}/results/qualitycontrol/fastqc/trimmed/{sample}_trimmed_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/trimmed/{sample}_trimmed_fastqc.zip"
    message: "Running FastQC on adapter-trimmed reads for sample {wildcards.sample} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/results/qualitycontrol/fastqc/trimmed/{sample}_trimmed_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.5.0/bio/fastqc"

# Rule: Run FastQC on quality-filtered reads
rule fastqc_quality_filtered:
    input:
        "{species}/processed/quality_filtered/{sample}_quality_filtered.fastq.gz"
    output:
        html="{species}/results/qualitycontrol/fastqc/quality_filtered/{sample}_quality_filtered_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/quality_filtered/{sample}_quality_filtered_fastqc.zip"
    message: "Running FastQC on quality-filtered reads for sample {wildcards.sample} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/results/qualitycontrol/fastqc/quality_filtered/{sample}_quality_filtered_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.5.0/bio/fastqc"

# Rule: Run FastQC on merged reads
rule fastqc_merged:
    input:
        merged="{species}/processed/merged/{individual}.fastq.gz"
    output:
        html="{species}/results/qualitycontrol/fastqc/merged/{individual}_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/merged/{individual}_fastqc.zip"
    message: "Running FastQC on merged reads for individual {wildcards.individual} in species {wildcards.species}"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/results/qualitycontrol/fastqc/merged/{individual}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.2.0/bio/fastqc"