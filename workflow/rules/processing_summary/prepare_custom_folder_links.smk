rule link_qualimap_for_multiqc:
    input:
        "{species}/results/{reference}/analytics/{individual}/qualimap"
    output:
        directory("{species}/results/summary/individual_level/{individual}/multiqc_custom_content/qualimap/{individual}_{reference}")
    shell:
        """
        mkdir -p $(dirname {output})
        cp -r {input} {output}
        """

rule copy_mapdamage_result_for_multiqc:
    input:
        GtoA3p  = "{species}/results/{reference}/analytics/{individual}/mapdamage/{individual}_{reference}.3pGtoA_freq.txt",
        CtoT5p  = "{species}/results/{reference}/analytics/{individual}/mapdamage/{individual}_{reference}.5pCtoT_freq.txt",
        lg_dist = "{species}/results/{reference}/analytics/{individual}/mapdamage/{individual}_{reference}.lgdistribution.txt",
    output:
        GtoA3p  = "{species}/results/summary/individual_level/{individual}/multiqc_custom_content/mapdamage/{individual}_{reference}/3pGtoA_freq.txt",
        CtoT5p  = "{species}/results/summary/individual_level/{individual}/multiqc_custom_content/mapdamage/{individual}_{reference}/5pCtoT_freq.txt",
        lg_dist = "{species}/results/summary/individual_level/{individual}/multiqc_custom_content/mapdamage/{individual}_{reference}/lgdistribution.txt",
    shell:
        """
        mkdir -p {wildcards.species}/results/summary/{wildcards.individual}/multiqc_custom_content/mapdamage/{wildcards.individual}_{wildcards.reference}
        cp {input.GtoA3p} {output.GtoA3p}
        cp {input.CtoT5p} {output.CtoT5p}
        cp {input.lg_dist} {output.lg_dist}
        """
