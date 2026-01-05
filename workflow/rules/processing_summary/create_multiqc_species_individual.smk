import textwrap

def create_multiqc_species_individual_input(wildcards):
    """Generate a list of input files for MultiQC report for a given species and its individuals."""

    species = wildcards.species
    individual = wildcards.individual

    #individuals = get_individuals_for_species(species)
    references = get_references_ids_for_species(species)

    file_list = []

    #for individual in individuals:
    # get samples for individual
    samples_of_individual = get_samples_for_species_individual(species, individual)

    for sample in samples_of_individual:

        # get raw read file paths
        raw_reads = get_raw_reads_for_sample(species, sample)

        # add fastp json reports
        # fastp trimming reports
        if len(raw_reads) == 2:
            file_list.append(f"{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.pe.json")
        else:
            file_list.append(f"{species}/results/reads/reads_trimmed/fastp_report/{sample}_trimmed.se.json")

        # fastp quality filtering reports
        file_list.append(f"{species}/results/reads/reads_quality_filtered/fastp_report/{sample}_quality_filtered.json")

        # add fastqc reports
        # raw reads fastqc
        #for read in raw_reads:
        #    file_list.append(f"{species}/results/reads/reads_raw/fastqc/{os.path.basename(read).replace('.fastq.gz','')}_raw_fastqc.zip")
        # trimmed reads fastqc
        #file_list.append(f"{species}/results/reads/reads_trimmed/fastqc/{sample}_trimmed_fastqc.zip")

        # eCMSD contamination analysis outputs
        file_list.append(f"{species}/results/summary/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_ReadLengths.png")
    
    file_list.append(f"{species}/results/contamination_analysis/ecmsd/{individual}_Mito_summary_genus_proportions_combined.tsv")
        
    # merged reads fastqc
    file_list.append(f"{species}/results/reads/reads_merged/fastqc/{individual}_merged_fastqc.zip")

    # bam analytics
    for reference in references:

        file_list.append(f"{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_reads_processing_summary.tsv",)
        file_list.append(f"{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_coverage_analysis.tsv")
        file_list.append(f"{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_depth_coverage_avg.csv")
        file_list.append(f"{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_coverage_summary.tsv")
        #file_list.append(f"{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_depth_coverage_max.csv")

        file_list.append(f"{species}/results/{reference}/plots/endogenous_reads/{species}_{reference}_raw_and_endogenous_reads_bar_chart.png")
        
        #file_list.append(f"{species}/results/{reference}/analytics/{individual}/picard_duplicates/{individual}_{reference}_metrics.txt")
        file_list.append(f"{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap")
        file_list.append(directory(f"{species}/results/{reference}/analytics/{individual}/qualimap"))
        file_list.append(f"{species}/results/{reference}/analytics/{individual}/samtools_stats/{individual}_{reference}_final.bam.stats")
        
        if config.get("pipeline", {}).get("reference_processing", {}).get("damage_rescaling", {}).get("execute", True) == True:
            file_list.append(directory(f"{species}/results/{reference}/analytics/{individual}/mapdamage/"))
        
        if config.get("pipeline", {}).get("reference_processing", {}).get("deduplication", {}).get("execute", True) == True:
            file_list.append(f"{species}/results/{reference}/analytics/{individual}/dedup/{individual}_{reference}_final.dedup.json")

    logger.debug(f"MultiQC input files for species {species}: {file_list}")

    return file_list

####################################################
# Snakemake rules
####################################################
rule create_multiqc_species_individual:
    input:
        create_multiqc_species_individual_input,
        config="{species}/results/summary/{individual}/{individual}_multiqc_config.yaml"
    output:
        "{species}/results/summary/{individual}_multiqc.html",
        directory("{species}/results/summary/{individual}/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
        use_input_files_only=True
    log:
        "{species}/results/summary/{individual}/multiqc.log",
    wrapper:
        "v7.9.0/bio/multiqc"

rule create_multiqc_species_individual_config:
    output:
        "{species}/results/summary/{individual}/{individual}_multiqc_config.yaml"
    script:
        "../../scripts/processing_summary/create_multiqc_species_individual_script_create_multiqc_species_individual_config.py"
        