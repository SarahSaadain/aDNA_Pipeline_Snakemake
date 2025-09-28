# Rule: Run eCMSD for contamination analysis
rule ecmsd_analysis:
    input:
        fastq = "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        outdir          = directory("{species}/results/contamination_analysis/ecmsd/{sample}"),
        summary         = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary.txt",
        paf             = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito.paf.gz",
        RMUS            = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary.RMUS.txt", 
        proportions     = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt",
        genus           = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary_genus.txt",
        readlength      = report("{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary_genus_ReadLengths.png",
                                category="contamination_analysis",
                                subcategory="ECMSD"),
        proportions_png = report("{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary_genus_Proportions.png",
                                category="contamination_analysis",
                                subcategory="ECMSD"),
    params:
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"]["executable"]
    threads: 15
    log:
        "{species}/logs/contamination_analysis/ecmsd/{sample}_ecmsd.log"
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"]["conda_env"]
    message: "Running eCMSD contamination analysis for {input.fastq}"
    shell:
        """
        bash {params.executable} \
            --fwd {input.fastq} \
            --out {output.outdir} \
            --threads {threads} \
            --prefix {wildcards.sample} \
            --Binsize 1000 \
            --RMUS-threshold 0.15 \
            --skip_environment  \
            --mapping_quality 20 \
            --taxonomic-hierarchy genus \
            --force \
        > {log} 2>&1
        """
