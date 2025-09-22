# 1 Map reads to reference (SAM output)
rule bwa_map_aDNA:
    input:
        reads="{species}/processed/merged/{individual}.fastq.gz",
        ref="{species}/raw/ref_genome/{ref_genome}.fa",
        index="{species}/raw/ref_genome/{ref_genome}.fa.bwt"
    output:
        sam=temp("{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.sam")
    message: "Mapping {input.reads} to reference genome {input.ref} for individual {wildcards.individual} in species {wildcards.species}"
    threads: workflow.cores 
    conda:
        "../../../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} > {output.sam}
        """

# 2 Convert SAM to BAM
rule sam_to_bam:
    input:
        "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.sam"
    output:
        bam=temp("{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.bam")
    message: "Converting SAM to BAM for {input}"
    threads: workflow.cores 
    wrapper:
        "v7.5.0/bio/samtools/view"

# 3 Sort BAM
rule sort_bam:
    input:
        "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.bam"
    output:
        "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam"
    message: "Sorting BAM file for {input}"
    log:
        "{species}/logs/{individual}_{ref_genome}_sort_bam.log",
    threads: workflow.cores 
    wrapper:
        "v7.5.0/bio/samtools/sort"

# 4 Index BAM
rule index_bam:
    input:
        "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam"
    output:
        "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam.bai"
    message: "Indexing BAM file for {input}"
    params:
        extra="",  # optional params string
    threads: workflow.cores
    wrapper:
        "v7.5.0/bio/samtools/index"
