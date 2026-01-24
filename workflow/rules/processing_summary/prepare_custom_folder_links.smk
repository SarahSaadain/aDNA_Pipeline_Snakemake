rule link_qualimap_for_multiqc:
    input:
        directory("{species}/results/{reference}/analytics/{individual}/qualimap")
    output:
        directory("{species}/results/summary/{individual}/multiqc_custom_content/qualimap/{individual}_{reference}")
    shell:
        """
        mkdir -p $(dirname {output})
        cp -r {input} {output}
        """

rule link_mapdamage_for_multiqc:
    input:
        directory("{species}/results/{reference}/analytics/{individual}/mapdamage"),
    output:
        directory("{species}/results/summary/{individual}/multiqc_custom_content/mapdamage/{individual}_{reference}/")
    shell:
        """
        mkdir -p $(dirname {output})
        cp -r {input} {output}
        """
