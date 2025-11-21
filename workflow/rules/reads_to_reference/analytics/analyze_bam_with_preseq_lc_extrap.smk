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

rule plot_analyze_bam_with_preseq_lc_extrap:
    input:
        "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap"
    output:
        "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap.pdf"
    log:
       "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap.pdf.log"
    shell:
        """
        python ../../../scripts/reads_to_reference/plotting/plot_preseq_complexity_curve.py {input} -o {output}
        """