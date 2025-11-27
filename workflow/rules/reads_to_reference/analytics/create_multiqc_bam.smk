####################################################
# Snakemake rules
####################################################

rule create_multiqc_bam:
    input:
        "{species}/results/{reference}/multiqc/"
    output:
        "{species}/results/{reference}/multiqc.html",
        directory("{species}/results/{reference}/analytics/multiqc_data")
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "{species}/results/{reference}/analytics/multiqc.log"
    wrapper:
        "v7.9.0/bio/multiqc"

rule create_multiqc_folder:
    input:
        insurance = "{species}/results/{reference}/analytics/request_multiqc_bam.ready"
    output:
        temp(directory("{species}/results/{reference}/multiqc/"))
    shell:
        "mkdir -p {output}"


rule validate_multiqc_preprocessing:
    input:
        picard_metrics = "{species}/results/{reference}/analytics/{individual}/picard_duplicates/{individual}_{reference}_metrics.txt",
        preseq_lc_extrap = "{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap",
        samtools_flagstat = "{species}/results/{reference}/analytics/{individual}/samtools_flagstats/{individual}_{reference}_final.bam.flagstats",
        samtools_stats = "{species}/results/{reference}/analytics/{individual}/samtools_stats/{individual}_{reference}_final.bam.stats",
        qualimap = "{species}/results/{reference}/analytics/{individual}/qualimap",
    output:
        preprocessing = temp("{species}/results/{reference}/analytics/{individual}/multiqc_preprocessing.ready")
    shell:
        "touch {output}"
    
rule request_multiqc_bam:
    input:
        preprocessing = lambda wildcards: expand(
            "{species}/results/{reference}/analytics/{individual}/multiqc_preprocessing.ready",
            species     = wildcards.species,
            reference   = wildcards.reference,
            individual  = get_individuals_for_species(wildcards.species),
        )
    output:
        ok_file = temp("{species}/results/{reference}/analytics/request_multiqc_bam.ready")
    shell:
        """
        echo 'Requesting preprocessing for multiqc report:'
        
        for f in {input.preprocessing}; do
            echo $f
        done

        touch {output.ok_file}
        """