# Rule: Map reads to reference genome using BWA
rule bwa_map_aDNA:
    # 1 Map reads to reference (SAM output)
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


# Rule: Convert SAM to BAM
rule sam_to_bam:
    # 2 Convert SAM to BAM
    input:
        "{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.sam"
    output:
        bam=temp("{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}.bam")
    message: "Converting SAM to BAM for {input}"
    threads: workflow.cores 
    wrapper:
        "v7.5.0/bio/samtools/view"


# Rule: Sort BAM file
rule sort_bam:
    # 3 Sort BAM
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


# Rule: Index BAM file
rule index_bam:
    # 4 Index BAM
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
