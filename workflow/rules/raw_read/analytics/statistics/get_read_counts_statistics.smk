import pandas as pd
import gzip

####################################################
# Python helper functions genereal
####################################################

def get_fastq_read_count(fastq_file):
    """
    counted reads in a FASTQ file (handles gzipped files).
    Each read is 4 lines, so counted lines and divide by 4.
    """

    logger.info(f"Counting reads in {fastq_file}")

    count = 0

    if fastq_file is None:
        return 0

    if fastq_file.endswith(".gz"):
        with gzip.open(fastq_file, "rt") as f:
            count = sum(1 for _ in f) // 4
    else:
        with open(fastq_file, "r") as f:
            count = sum(1 for _ in f) // 4

    logger.info(f"Found {count} reads in {fastq_file}")
    return count

####################################################
# Snakemake rules
####################################################

# Rule: Count reads in raw FASTQ files
rule count_reads_raw:
    input:
        fastq=lambda wc: get_raw_reads_for_sample(wc.species, wc.sample), 
    output:
        counted="{species}/processed/reads/statistics/{sample}_raw.count"
    message: "Counting reads in raw FASTQ file(s) {input.fastq}"
    conda:
        "../../../../envs/python_and_r.yaml",
    run:
        files = input.fastq if isinstance(input.fastq, list) else [input.fastq]
        count = sum(get_fastq_read_count(f) for f in files)
        with open(output.counted, "w") as f:
            f.write(str(count))

# Rule: Count reads in trimmed FASTQ files
rule count_reads_trimmed:
    input:
        fastq="{species}/processed/reads/reads_trimmed/{sample}_trimmed_final.fastq.gz"
    output:
        counted="{species}/processed/reads/statistics/{sample}_trimmed.count"
    message: "Counting reads in {input.fastq}"
    conda:
        "../../../../envs/python_and_r.yaml",
    run:
        count = get_fastq_read_count(input.fastq)
        with open(output.counted, "w") as f:
            f.write(str(count))

# Rule: Count reads in quality-filtered FASTQ files
rule count_reads_quality_filtered:
    input:
        fastq="{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz"
    output:
        counted="{species}/processed/reads/statistics/{sample}_quality_filtered.count"
    message: "Counting reads in {input.fastq}"
    conda:
        "../../../../envs/python_and_r.yaml",
    run:
        count = get_fastq_read_count(input.fastq)
        with open(output.counted, "w") as f:
            f.write(str(count))

# Rule: Combine read counts per sample
rule combine_counts_per_sample:
    input:
        raw_reads="{species}/processed/reads/statistics/{sample}_raw.count",
        trimmed_reads="{species}/processed/reads/statistics/{sample}_trimmed.count",
        quality_filtered_reads="{species}/processed/reads/statistics/{sample}_quality_filtered.count"
    output:
        counts="{species}/processed/reads/statistics/{sample}_reads_counts.csv"
    message: "Combining read counts for sample {wildcards.sample}"
    conda:
        "../../../../envs/python_and_r.yaml",
    run:

        #print(input.raw_reads)
        #print(input.trimmed_reads)
        #print(input.quality_filtered_reads)
        #print(output.counts)
        #print(wildcards.sample)

        with open(input.raw_reads, "r") as f:
            raw = int(f.read())
        with open(input.trimmed_reads, "r") as f:
            trimmed = int(f.read())
        with open(input.quality_filtered_reads, "r") as f:
            quality_filtered = int(f.read())

        parts = wildcards.sample.split("_")
        individual = parts[0] if len(parts) > 0 else "N/A"      

        with open(output.counts, "w") as f:
            f.write(f"{wildcards.sample},{individual},{str(raw)},{str(trimmed)},{str(quality_filtered)}\n")

# Rule: Combine read counts per species
rule combine_counts_per_species:
    input:
        lambda wildcards: expand("{species}/processed/reads/statistics/{sample}_reads_counts.csv",
            sample=get_sample_ids_for_species(wildcards.species),
            species=wildcards.species)
    output:
        counts="{species}/results/reads/statistics/{species}_reads_counts.csv"
    conda:
        "../../../../envs/python_and_r.yaml",
    run:
        data = []
        
        pd.concat([pd.read_csv(f, header=None, names=["reads_file", "individual", "raw_count", "adapter_removed_count", "quality_filtered_count"]) for f in input], ignore_index=True).to_csv(output.counts, index=False)