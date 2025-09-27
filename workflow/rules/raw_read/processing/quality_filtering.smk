def get_quality_filtered_input_read(wildcards):
    # Determine if the sample is paired-end or single-end
    reads = get_adapter_removal_input_reads(wildcards)
    if len(reads) == 2:
        # Paired-end: use the merged reads from fastp_pe
        return [f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.pe.fastq.gz"]
    else:
        # Single-end: use the trimmed reads from fastp_se
        return [f"{wildcards.species}/processed/reads/reads_trimmed/{wildcards.sample}_trimmed.se.fastq.gz"]
 
# Rule: Quality filtering of reads using fastp
rule quality_filter:
    input:
        sample=get_quality_filtered_input_read
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
        "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.log",
    params:
        extra="--disable_adapter_trimming --qualified_quality_phred 15 --length_required 15 --unqualified_percent_limit 40 --n_base_limit 5"
    threads: workflow.cores
    wrapper:
        "v7.5.0/bio/fastp"