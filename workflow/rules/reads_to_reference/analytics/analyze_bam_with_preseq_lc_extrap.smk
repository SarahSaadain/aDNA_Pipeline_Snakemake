# NOTE: preseq must be run on the raw mapped BAM before any duplicate removal or damage rescaling.
# Running preseq on a deduplicated and/or mapDamage-rescaled BAM removes duplication information
# and violates preseq assumptions, leading to errors such as:
# “ERROR: Saturation expected at double initial sample size. Unable to extrapolate”.
# Library complexity estimates from deduplicated or rescaled BAMs are invalid and cannot be recovered.

rule analyze_bam_with_preseq_lc_extrap:
    input:
        # bam files containing duplicates and sorted by using bamtools or samtools sort.
        "{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam"
    output:
        "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap"
    params:
        "-v"   #optional parameters
    log:
       "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap.log"
    wrapper:
        "v2.10.0/bio/preseq/lc_extrap"