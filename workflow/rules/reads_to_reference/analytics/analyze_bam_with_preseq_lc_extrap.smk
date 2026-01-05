rule analyze_bam_with_preseq_lc_extrap:
    input:
        # bam files containing duplicates and sorted by using bamtools or samtools sort.
        "{species}/processed/{reference}/mapped/{individual}_{reference}_final.bam"
    output:
        "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap"
    params:
        "-v"   #optional parameters
    log:
       "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap.log"
    wrapper:
        "v2.10.0/bio/preseq/lc_extrap"