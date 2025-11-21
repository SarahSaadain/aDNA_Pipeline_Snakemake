####################################################
# Snakemake rules
####################################################

rule samtools_flagstat:
    input:
        "{species}/processed/{reference}/mapped/{individual}_{reference}_final.bam"
    output:
        "{species}/results/{reference}/analytics/{individual}/samtools_flagstats/{individual}_{reference}_final.bam.flagstats"
    log:
        "{species}/results/{reference}/analytics/{individual}/samtools_flagstats/{individual}_{reference}_final.bam.flagstats.log"
    params:
        extra="",  # optional params string
    wrapper:
        "v7.6.0/bio/samtools/flagstat"