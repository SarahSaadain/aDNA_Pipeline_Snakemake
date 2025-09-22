rule ecmsd_analysis:
    input:
        fastq = "{species}/processed/quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        outdir = directory("{species}/results/qualitycontrol/ecmsd/{sample}"),
        logdir = directory("{species}/results/qualitycontrol/ecmsd/{sample}/logs"),
        summary = "{species}/results/qualitycontrol/ecmsd/{sample}/mapping/Mito_summary.txt",
        mito_paf = "{species}/results/qualitycontrol/ecmsd/{sample}/mapping/Mito.paf.gz"
    params:
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"]["executable"]
    threads: workflow.cores
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"]["conda_env"]
    message: "Running eCMSD contamination analysis for {input.fastq}"
    shell:
        """
        bash {params.executable} \
            --fwd {input.fastq} \
            --out {output.outdir} \
            --threads {threads} \
            --Binsize 1000 \
            --RMUS-threshold 0.15 \
            --mapping_quality 20 \
            --taxonomic-hierarchy genus \
            --force
        """
