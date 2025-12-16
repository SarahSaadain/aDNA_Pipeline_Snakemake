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
        #file_list.append(f"{species}/results/reads/reads_quality_filtered/fastp_report/{sample}_quality_filtered.json")

        # add fastqc reports
        # raw reads fastqc
        for read in raw_reads:
            file_list.append(f"{species}/results/reads/reads_raw/fastqc/{os.path.basename(read).replace('.fastq.gz','')}_raw_fastqc.zip")
        # trimmed reads fastqc
        file_list.append(f"{species}/results/reads/reads_trimmed/fastqc/{sample}_trimmed_fastqc.zip")
        # quality filtered reads fastqc
        #file_list.append(f"{species}/results/reads/reads_quality_filtered/fastqc/{sample}_quality_filtered_fastqc.zip")

        # eCMSD contamination analysis outputs
        file_list.append(f"{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_ReadLengths.png")
        file_list.append(f"{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_Proportions.png")
        #file_list.append(f"{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary.txt")
        file_list.append(f"{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_proportions.txt")
        

        # merged reads fastqc
        file_list.append(f"{species}/results/reads/reads_merged/fastqc/{individual}_merged_fastqc.zip")
    
        # bam analytics
        for reference in references:

            file_list.append(f"{species}/results/{reference}/analytics/{individual}/endogenous/{individual}_{reference}.endogenous.csv")
            file_list.append(f"{species}/results/{reference}/analytics/{individual}/coverage/{individual}_{reference}_coverage_analysis.csv")

            file_list.append(f"{species}/results/{reference}/plots/endogenous_reads/{species}_{reference}_raw_and_endogenous_reads_bar_chart.png")
            
            file_list.append(f"{species}/results/{reference}/analytics/{individual}/picard_duplicates/{individual}_{reference}_metrics.txt")
            file_list.append(f"{species}/results/{reference}/analytics/{individual}/preseq/{individual}_{reference}.lc_extrap")
            file_list.append(directory(f"{species}/results/{reference}/analytics/{individual}/qualimap"))
            #file_list.append(f"{species}/results/{reference}/analytics/{individual}/samtools_flagstats/{individual}_{reference}_final.bam.flagstats")
            file_list.append(f"{species}/results/{reference}/analytics/{individual}/samtools_stats/{individual}_{reference}_final.bam.stats")
            
            if config.get("pipeline", {}).get("reference_processing", {}).get("damage_rescaling", {}).get("execute", True) == True:
                file_list.append(directory(f"{species}/results/{reference}/analytics/{individual}/mapdamage/"))
            
            if config.get("pipeline", {}).get("reference_processing", {}).get("deduplication", {}).get("execute", True) == True:
                file_list.append(f"{species}/results/{reference}/analytics/{individual}/dedup/{individual}_{reference}_sorted.dedup.json")

    logger.debug(f"MultiQC input files for species {species}: {file_list}")

    return file_list

def generate_multiqc_config_for_species_and_individual(species, individual):

    return (
        f"# MultiQC configuration for species: {species}, individual: {individual}\n"
        f"title: 'Snakepipe Report for species {species}, individual {individual}'\n"
       # f"subtitle: \"Results overview for species {species}, individual {individual} from SnakePipe aDNA Pipeline\"\n"
        f"intro_text: 'This MultiQC report summarises the results for species {species}, individual {individual} processed with the SnakePipe aDNA Pipeline.'\n"
        "\n"
        "custom_data:\n"
        #"   ecmsd:\n"
        #"       section_name: 'ECMSD Contamination Analysis'\n"
        #"       parent_id: ecmsd\n"
        "   ecmsd_proportions:\n"
        "       parent_id: ecmsd\n"
        "       parent_name: 'ECMSD Contamination Analysis'\n"
        "       section_name: 'ECMSD Proportions Plot'\n"
        "       description: 'Relative Abundance of Mitochondrial Reads Assigned to Different Genera.'\n"
        "       file_format: 'png'\n"
         "\n"
        "   ecmsd_proportions_table:\n"
        "       parent_id: ecmsd\n"
        "       parent_name: 'ECMSD Contamination Analysis'\n"
        "       section_name: 'ECMSD Proportions Table'\n"
        "       description: 'This table visualizes the proportions of reads mapped to different Mitochondrial genomes to check for contamination.'\n"
        "       file_format: 'tsv'\n"
        "       plot_type: 'table'\n"
        "\n"
        "   ecmsd_readlengths:\n"
        "       parent_id: ecmsd\n"
        "       parent_name: 'ECMSD Contamination Analysis'\n"
        "       section_name: 'ECMSD Read Lengths'\n"
        "       description: 'Distribution of read lengths across genera to evaluate contamination and distinguish ancient from modern read.'\n"
        "       file_format: 'png'\n"
        "\n"
        "   depth_table:\n"        
        "       parent_id: depth\n"
        "       parent_name: 'Depth of Coverage'\n"
        "       section_name: 'Depth of Coverage'\n"
        "       description: 'Depth of coverage across the reference genome.'\n"
        "       file_format: 'csv'\n"
         "\n"
        "   endogenous_table:\n"
        "       parent_id: endogenous\n"
        "       parent_name: 'Endogenous Reads'\n"
        "       section_name: 'Endogenous Reads Table'\n"
        "       description: 'This table visualizes the proportions of reads mapped to the reference.'\n"
        "       file_format: 'csv'\n"
        "\n"
        "   endogenous_plot:\n"
        "       parent_id: endogenous\n"
        "       parent_name: 'Endogenous Reads'\n"
        "       section_name: 'Endogenous Reads Bar Chart'\n"
        "       description: 'Bar chart showing the number of endogenous reads.'\n"
        "       file_format: 'png'\n"
        "\n"
        "sp:\n"
        "   ecmsd_readlengths:\n"
        "       fn: \"*_Mito_summary_genus_ReadLengths.png\"\n"
        "   ecmsd_proportions:\n"
        "       fn: \"*_Mito_summary_genus_Proportions.png\"\n"
        "   ecmsd_proportions_table:\n"
        "       fn: \"*_Mito_summary_genus_proportions.txt\"\n"
        "   depth_table:\n"
        "       fn: '*_coverage_analysis.csv'\n"
        "   endogenous_table:\n"
        "       fn: '*.endogenous.csv'\n"
        "   endogenous_plot:\n"
        "       fn: '*_raw_and_endogenous_reads_bar_chart.png'\n"
        "\n"
        "ignore_images: false\n"
        "module_order:\n"
        "   - fastqc:\n"
        "       name: 'FastQC (raw)'\n"
        "       anchor: 'fastqc_raw'\n"
        "       info: 'This section of the report shows FastQC results for the raw reads.'\n"
        "       path_filters:\n"
        "         - '*_raw_fastqc.zip'\n"
        "   - fastqc:\n"
        "       name: 'FastQC (adapter removed)'\n"
        "       anchor: 'fastqc_trimmed'\n"
        "       info: 'This section of the report shows FastQC results after adapter trimming.'\n"
        "       target: ''\n"
        "       path_filters:\n"
        "         - '*_trimmed_fastqc.zip'\n"
        "   - fastqc:\n"
        "       name: 'FastQC (quality filtered)'\n"
        "       anchor: 'fastqc_filtered'\n"
        "       info: 'This section of the report shows FastQC results after quality filtering.'\n"
        "       target: ''\n"
        "       path_filters:\n"
        "         - '*_quality_filtered_fastqc.zip'\n"
        "   - fastqc:\n"
        "       name: 'FastQC (merged)'\n"
        "       anchor: 'fastqc_merged'\n"
        "       info: 'This section of the report shows FastQC results after merging.'\n"
        "       target: ''\n"
        "       path_filters:\n"
        "         - '*_merged_fastqc.zip'\n"
        "   - fastp:\n"
        "       name: 'Fastp (adapter removal)'\n"
        "       anchor: 'fastp_adapter_removal'\n"
        "       info: 'This section of the report shows FastP results of adapter trimming.'\n"
        "       target: ''\n"
        "       path_filters:\n"
        "         - '*_trimmed.*.json'\n"
        "   - fastp:\n"
        "       name: 'Fastp (quality filtering)'\n"
        "       anchor: 'fastp_quality_filtered'\n"
        "       info: 'This section of the report shows FastP results of quality filtering.'\n"
        "       target: ''\n"
        "       path_filters:\n"
        "         - '*_quality_filtered.json'\n"
        "report_section_order:\n"
        "   fastqc_raw:\n"
        "       order: 130\n"
        "   fastp_adapter_removal:\n"
        "       order: 120\n"
        "   fastqc_trimmed:\n"
        "       order: 110\n"
        "   fastp_quality_filtered:\n"
        "       order: 100\n"
        "   fastqc_filtered:\n"
        "       order: 90\n"
        "   fastqc_merged:\n"
        "       order: 80\n"
        "   ecmsd:\n"
        "       order: 70\n"
        "   endogenous:\n"
        "       order: 62\n"
        "   depth:\n"
        "       order: 61\n"
        "   qualimap:\n"
        "       order: 60\n"
        "   samtools:\n"
        "       order: 50\n"
        "   dedup:\n"
        "       order: 40\n"
        "   picard:\n"
        "       order: 30\n"
        "   preseq:\n"
        "       order: 20\n"
        "   mapdamage:\n"
        "       order: 10\n"
    )

####################################################
# Snakemake rules
####################################################
rule create_multiqc_species_individual:
    input:
        create_multiqc_species_individual_input,
        config="{species}/results/analytics/{individual}/{individual}_multiqc_config.yaml"
    output:
        "{species}/results/analytics/{individual}_multiqc.html",
        directory("{species}/results/analytics/{individual}/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
        use_input_files_only=True
    log:
        "{species}/results/analytics/{individual}/multiqc.log",
    wrapper:
        "v7.9.0/bio/multiqc"

rule create_multiqc_species_individual_config:
    output:
        "{species}/results/analytics/{individual}/{individual}_multiqc_config.yaml"
    run:
        species = wildcards.species
        individual = wildcards.individual

        # Generate the MultiQC config content
        config_content = generate_multiqc_config_for_species_and_individual(species, individual)

        # Write to output file
        with open(output[0], 'w') as f:
            f.write(config_content)

rule ecmsd_for_multiqc_report:
    input:
        summary = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary.txt",
        length = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary_genus_ReadLengths.png",
        proportions = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary_genus_Proportions.png",
        proportions_txt = "{species}/results/contamination_analysis/ecmsd/{sample}/mapping/{sample}_Mito_summary_genus_proportions.txt"
    output:
        summary = "{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary.txt",
        length = "{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_ReadLengths.png",
        proportions = "{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_Proportions.png",
        proportions_txt = "{species}/results/analytics/{individual}/multiqc_custom_content/{sample}/{sample}_Mito_summary_genus_proportions.txt"
    shell:
        """
        cp {input.summary} {output.summary}
        cp {input.length} {output.length}
        cp {input.proportions} {output.proportions}
        cp {input.proportions_txt} {output.proportions_txt}
        """