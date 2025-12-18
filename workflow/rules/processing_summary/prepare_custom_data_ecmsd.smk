rule ecmsd_for_multiqc_report:
    input:
        summary = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary.txt",
        length = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_ReadLengths.png",
        #proportions = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_Proportions.png",
        #proportions_txt = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt"
    output:
        summary = "{species}/results/summary/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary.txt",
        length = "{species}/results/summary/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_ReadLengths.png",
        #proportions = "{species}/results/summary/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_Proportions.png",
        #proportions_txt = "{species}/results/summary/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_proportions.txt"
    shell:
        """
        cp {input.summary} {output.summary}
        cp {input.length} {output.length}
        """