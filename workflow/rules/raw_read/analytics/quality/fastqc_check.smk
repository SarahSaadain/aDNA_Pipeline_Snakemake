rule fastqc_raw:
    input:
        "{species}/raw/reads/{sample}.fastq.gz"
    output:
        html="{species}/results/qualitycontrol/fastqc/raw/{sample}_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/raw/{sample}_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/results/qualitycontrol/fastqc/raw/{sample}_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.2.0/bio/fastqc"

rule fastqc_adapter_removed:
    input:
        "{species}/processed/trimmed/{sample}_trimmed.fastq.gz"
    output:
        html="{species}/results/qualitycontrol/fastqc/trimmed/{sample}_trimmed_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/trimmed/{sample}_trimmed_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/results/qualitycontrol/fastqc/trimmed/{sample}_trimmed_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.2.0/bio/fastqc"

rule fastqc_quality_filtered:
    input:
        "{species}/processed/quality_filtered/{sample}_quality_filtered.fastq.gz"
    output:
        html="{species}/results/qualitycontrol/fastqc/quality_filtered/{sample}_quality_filtered_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/quality_filtered/{sample}_quality_filtered_fastqc.zip"
    params:
        extra="--quiet",
        mem_overhead_factor=0.1,
    log:
        "{species}/results/qualitycontrol/fastqc/quality_filtered/{sample}_quality_filtered_fastqc.log"
    threads: 1
    resources:
        mem_mb = 1024,
    wrapper:
        "v7.2.0/bio/fastqc"

rule fastqc_merged:
    input:
        merged="{species}/processed/merged/{individual}.fastq.gz"
    output:
        html="{species}/results/qualitycontrol/fastqc/merged/{individual}_fastqc.html",
        zip="{species}/results/qualitycontrol/fastqc/merged/{individual}_fastqc.zip"
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