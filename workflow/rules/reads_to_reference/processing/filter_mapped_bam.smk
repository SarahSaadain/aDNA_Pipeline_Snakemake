####################################################
# Snakemake rules
####################################################

rule deduplicate_bam_with_dedup:
    input:
        bam="{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam"
    output:
        dedup_folder=directory("{species}/processed/{reference}/mapped/deduplication_{individual}/"),
        dedup_bam="{species}/processed/{reference}/mapped/deduplication_{individual}/{individual}_{reference}_sorted_rmdup.bam",
        dedup_hist="{species}/processed/{reference}/mapped/deduplication_{individual}/{individual}_{reference}_sorted.hist",
        dedup_json="{species}/processed/{reference}/mapped/deduplication_{individual}/{individual}_{reference}_sorted.dedup.json",
    message:
        "Deduplicating BAM file for {input.bam} using dedup for individual {wildcards.individual} in species {wildcards.species}",
    conda:
        "../../../envs/dedup.yaml"
    shell:
        """
        mkdir -p {output.dedup_folder}
        dedup --input {input.bam} --merged --output {output.dedup_folder}
        """