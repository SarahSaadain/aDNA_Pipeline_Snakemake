####################################################
# Snakemake rules
####################################################

# Rule: Map reads to reference using BWA
rule map_reads_to_reference:
    input:
        reads=["{species}/processed/reads/reads_merged/{individual}.fastq.gz"],
        idx=multiext("{species}/raw/ref/{reference}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        temp("{species}/processed/{reference}/mapped/{individual}_{reference}_unsorted.bam")
    log:
        "{species}/logs/{reference}/mapped/{individual}_{reference}.bam.log"
    threads: 15
    wrapper:
        "v7.6.0/bio/bwa/mem"

# Rule: Sort BAM file
rule sort_mapped_reads_bam:
    # 3 Sort BAM
    input:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_unsorted.bam"
    output:
        temp("{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam")
    message: "Sorting BAM file for {input}"
    log:
        "{species}/logs/{reference}/mapped/{individual}_{reference}_sorted_bam.log",
    threads: 10
    wrapper:
        "v7.5.0/bio/samtools/sort"

# Rule: Index BAM file
rule index_mapped_sorted_reads_bam:
    # 4 Index BAM
    input:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam" 
    output:
        temp("{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam.bai")
    message: "Indexing BAM file for {input}"
    params:
        extra="",  # optional params string
    threads: 5
    wrapper:
        "v7.5.0/bio/samtools/index"
