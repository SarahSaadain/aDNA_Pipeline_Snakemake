####################################################
# Snakemake rules
####################################################

# Rule: Run eCMSD for contamination analysis
rule ecmsd_analyze_contamination:
    input:
        fastq = "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        outdir          = directory("{species}/results/contamination_analysis/ecmsd/{individual}/{sample}"),
        summary         = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary.txt",
        paf             = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito.paf.gz",
        RMUS            = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary.RMUS.txt", 
        proportions     = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt",
        genus           = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus.txt",
        readlength      = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_ReadLengths.png",
        proportions_png = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_Proportions.png",
    params:
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"]["executable"]
    threads: 15
    log:
        "{species}/logs/contamination_analysis/ecmsd/{individual}/{sample}/{sample}_ecmsd.log"
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"].get("conda_env", "../../../../envs/ecmsd.yaml")
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

rule ecmsd_merge_genus_per_individual:
    input:
        lambda wildcards: expand("{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt",
            sample=get_samples_for_species_individual(wildcards.species, wildcards.individual),
            species=wildcards.species,
            individual=wildcards.individual
            )
    output:
        "{species}/results/contamination_analysis/ecmsd/{individual}_Mito_summary_genus_proportions_combined.tsv"
    script:
        "../../../../scripts/raw_reads/analytics/contamination/check_contamination_ecmsd_script_ecmsd_merge_genus_per_individual.py"
