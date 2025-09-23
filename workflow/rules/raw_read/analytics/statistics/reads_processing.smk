import pandas as pd
import gzip
#import input_manager as im

def get_fastq_read_count(fastq_file):
    """
    counted reads in a FASTQ file (handles gzipped files).
    Each read is 4 lines, so counted lines and divide by 4.
    """
    if fastq_file is None:
        return 0
    if fastq_file.endswith(".gz"):
        with gzip.open(fastq_file, "rt") as f:
            return sum(1 for _ in f) // 4
    else:
        with open(fastq_file, "r") as f:
            return sum(1 for _ in f) // 4

#def get_samples_per_species(species):
#    im.get_input_reads_processing_raw(species)

# Rule: Count reads in raw FASTQ files
rule count_reads_raw:
    input:
        fastq=lambda wildcards: get_files_in_folder_matching_pattern(os.path.join(wildcards.species, "raw", "reads"), f"{wildcards.sample}*R1*.fastq.gz")[0]
    output:
        counted=temp("{species}/processed/qualitycontrol/statistics/{sample}_raw_reads.count")
    message: "Counting reads in raw FASTQ file {input.fastq}"
    run:
        count = get_fastq_read_count(input.fastq)
        with open(output.counted, "w") as f:
            f.write(str(count) + "\n")

# Rule: Count reads in trimmed FASTQ files
rule count_reads_trimmed:
    input:
        fastq="{species}/processed/trimmed/{sample}_trimmed.fastq.gz"
    output:
        counted=temp("{species}/processed/qualitycontrol/statistics/{sample}_trimmed_reads.count")
    message: "Counting reads in {input.fastq}"
    run:
        count = get_fastq_read_count(input.fastq)
        with open(output.counted, "w") as f:
            f.write(str(count) + "\n")

# Rule: Count reads in quality-filtered FASTQ files
rule count_reads_quality_filtered:
    input:
        fastq="{species}/processed/quality_filtered/{sample}_quality_filtered.fastq.gz"
    output:
        counted=temp("{species}/processed/qualitycontrol/statistics/{sample}_quality_filtered_reads.count")
    message: "Counting reads in {input.fastq}"
    run:
        count = get_fastq_read_count(input.fastq)
        with open(output.counted, "w") as f:
            f.write(str(count) + "\n")

# Rule: Combine read counts per sample
rule combine_counts_per_sample:
    input:
        raw_reads="{species}/processed/qualitycontrol/statistics/{sample}_raw_reads.count",
        trimmed_reads="{species}/processed/qualitycontrol/statistics/{sample}_trimmed_reads.count",
        quality_filtered_reads="{species}/processed/qualitycontrol/statistics/{sample}_quality_filtered_reads.count"
    output:
        counts=temp("{species}/processed/qualitycontrol/statistics/{sample}_reads_processing.csv")
    message: "Combining read counts for sample {wildcards.sample}"
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
            f.write(wildcards.sample + "," + individual + "," + str(raw) + "," + str(trimmed) + "," + str(quality_filtered) + "\n")

# Rule: Combine read counts per species
rule combine_counts_per_species:
    input:
        lambda wildcards: expand("{species}/processed/qualitycontrol/statistics/{sample}_reads_processing.csv",
            sample=get_sample_ids_for_species(wildcards.species),
            species=wildcards.species)
    output:
        counts="{species}/results/qualitycontrol/statistics/{species}_reads_processing.csv"
    run:
        data = []
        
        pd.concat([pd.read_csv(f, header=None, names=["reads_file", "individual", "raw_count", "adapter_removed_count", "quality_filtered_count"]) for f in input], ignore_index=True).to_csv(output.counts, index=False)