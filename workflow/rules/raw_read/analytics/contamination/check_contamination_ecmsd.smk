####################################################
# Snakemake rules
####################################################

# Rule: Run eCMSD for contamination analysis
rule ecmsd_analyze_contamination:
    input:
        fastq = "{species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz",
    output:
        #outdir          = directory("{species}/results/contamination_analysis/ecmsd/{individual}/{sample}"),
        summary         = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary.txt",
        paf             = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito.paf.gz",
        RMUS            = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary.RMUS.txt", 
        proportions     = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt",
        genus           = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus.txt",
        readlength      = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_ReadLengths.png",
        proportions_png = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_Proportions.png",
    params:
        executable = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"]["executable"],
        binsize = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"].get("Binsize", 1000),
        rmus_threshold = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"].get("RMUS_threshold", 0.15),
        mapping_quality = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"].get("mapping_quality", 20),
        taxonomic_hierarchy = config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"].get("taxonomic_hierarchy", "genus"),
    threads: 15
    log:
        "{species}/logs/contamination_analysis/ecmsd/{individual}/{sample}/{sample}_ecmsd.log"
    conda:
        config["pipeline"]["raw_reads_processing"]["contamination_analysis"]["tools"]["ecmsd"]["settings"].get("conda_env", "../../../../envs/ecmsd.yaml")
    message: "Running eCMSD contamination analysis for {input.fastq}"
    shell:
        """
        # get 2 folders up for output from summary file
        outdir=$(dirname $(dirname {output.summary}))

        # create output folder
        mkdir -p "$outdir"

        echo "Running eCMSD for sample {wildcards.sample}"
        echo "Input FASTQ: {input.fastq}"
        echo "Output folder: $outdir"

        bash {params.executable} \
            --fwd {input.fastq} \
            --out "$outdir" \
            --threads {threads} \
            --prefix {wildcards.sample} \
            --Binsize {params.binsize} \
            --RMUS-threshold {params.rmus_threshold} \
            --skip_environment \
            --mapping_quality {params.mapping_quality} \
            --taxonomic-hierarchy {params.taxonomic_hierarchy} \
            --force \
            > {log} 2>&1
        """

rule ecmsd_merge_hits_per_individual:
    input:
        lambda wildcards: expand("{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt",
            sample=get_samples_for_species_individual(wildcards.species, wildcards.individual),
            species=wildcards.species,
            individual=wildcards.individual
            )
    output:
        "{species}/results/contamination_analysis/ecmsd/{individual}_Mito_summary_genus_hits_combined.tsv"
    script:
        "../../../../scripts/raw_reads/analytics/contamination/check_contamination_ecmsd_script_ecmsd_merge_hits_per_individual.py"

rule ecmsd_analyze_proportions:
    input:
        report = "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt",
        count_reads = "{species}/processed/reads/statistics/{sample}_quality_filtered.count"
    output:
        "{species}/results/contamination_analysis/ecmsd/{individual}/{sample}/pipeline/{sample}_ecmsd_proportions.tsv"
    params:
        sample = "{sample}"
    script:
        "../../../../scripts/raw_reads/analytics/contamination/check_contamination_ecmsd_script_ecmsd_analyze_proportions.py"
