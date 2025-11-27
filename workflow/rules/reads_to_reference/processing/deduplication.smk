####################################################
# Snakemake rules
####################################################

rule deduplicate_bam_with_dedup:
    input:
        bam="{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam"
    output:
        dedup_folder= directory("{species}/processed/{reference}/deduplication/{individual}/"),
        dedup_bam   = "{species}/processed/{reference}/deduplication/{individual}/{individual}_{reference}_sorted_rmdup.bam",
        dedup_hist  = "{species}/processed/{reference}/deduplication/{individual}/{individual}_{reference}_sorted.hist",
        dedup_json  = "{species}/processed/{reference}/deduplication/{individual}/{individual}_{reference}_sorted.dedup.json",
    message:
        "Deduplicating BAM file for {input.bam} using dedup for individual {wildcards.individual} in species {wildcards.species}",
    log: 
        "{species}/processed/{reference}/deduplication/{individual}_dedup.log",
    resources:
        mem_mb = 20000
    conda:
        "../../../envs/dedup.yaml"
    shell:
        """
        mkdir -p {output.dedup_folder}
        dedup -Xms5g -Xmx20g --input {input.bam} --merged --output {output.dedup_folder} > {log}
        """

rule move_deduplicated_to_mapped:
    input:
        dedup_folder =  "{species}/processed/{reference}/deduplication/{individual}/",
        dedup_bam = "{species}/processed/{reference}/deduplication/{individual}/{individual}_{reference}_sorted_rmdup.bam",
        dedup_hist  = "{species}/processed/{reference}/deduplication/{individual}/{individual}_{reference}_sorted.hist",
        dedup_json  = "{species}/processed/{reference}/deduplication/{individual}/{individual}_{reference}_sorted.dedup.json",
    output:
        moved_bam = temp("{species}/processed/{reference}/mapped/{individual}_{reference}_unsorted_dedupped.bam"),
        moved_folder = directory("{species}/results/{reference}/analytics/{individual}/dedup/"),
        moved_hist = "{species}/results/{reference}/analytics/{individual}/dedup/{individual}_{reference}_sorted.hist",
        moved_json = "{species}/results/{reference}/analytics/{individual}/dedup/{individual}_{reference}_sorteed.dedup.json",
    shell:
        """
        mv {input.dedup_bam} {output.moved_bam}
        mkdir -p {output.moved_folder}
        mv {input.dedup_hist} {output.moved_hist}
        mv {input.dedup_json} {output.moved_json}

        rm -r {input.dedup_folder}
        """

# Rule: Sort BAM file
rule sort_mapped_dedupped_reads_bam:
    # 3 Sort BAM
    input:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_unsorted_dedupped.bam"
    output:
        temp("{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped.bam")
    message: "Sorting BAM file for {input}"
    log:
        "{species}/logs/{reference}/mapped/{individual}_{reference}_sorted_bam.log",
    threads: 10
    wrapper:
        "v7.5.0/bio/samtools/sort"


# Rule: Index BAM file
rule index_mapped_sorted_dedupped_reads_bam:
    # 4 Index BAM
    input:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped.bam" 
    output:
        temp("{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped.bam.bai")
    message: "Indexing BAM file for {input}"
    params:
        extra="",  # optional params string
    threads: 5
    wrapper:
        "v7.5.0/bio/samtools/index"