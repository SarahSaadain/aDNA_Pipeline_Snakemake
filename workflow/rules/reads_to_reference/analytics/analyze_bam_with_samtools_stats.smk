####################################################
# Snakemake rules
####################################################

# Rule: Generate samtools stats for BAM files
rule analyze_bam_with_samtools_stats:
    input:
        bam="{species}/processed/{reference}/mapped/{individual}_{reference}_final.bam"
    output:
        "{species}/results/{reference}/analytics/{individual}/samtools_stats/{individual}_{reference}_final.bam.stats"
    message: "Generating samtools stats for {input.bam}"
    log:
        "{species}/logs/{reference}/statistics/{individual}/{individual}_{reference}_final.bam.stats.log"
    wrapper:
        "v7.2.0/bio/samtools/stats"