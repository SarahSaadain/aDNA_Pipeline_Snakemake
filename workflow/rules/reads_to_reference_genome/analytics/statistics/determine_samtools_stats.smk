# Rule: Generate samtools stats for BAM files
rule samtools_stats:
    input:
        bam="{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.rescaled.bam"
    output:
        "{species}/results/{ref_genome}/statistics/{individual}/{individual}_{ref_genome}.bam.stats"
    message: "Generating samtools stats for {input.bam}"
    log:
        "{species}/results/{ref_genome}/statistics/{individual}/{individual}_{ref_genome}.bam.stats.log"
    wrapper:
        "v7.2.0/bio/samtools/stats"