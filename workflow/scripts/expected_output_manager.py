# =================================================================================================
#     Input Manager Utility Functions for aDNA Pipeline
# =================================================================================================
# This script provides helper functions for managing input files and sample metadata in the pipeline.
# Functions include file discovery, reference handling, and sample identification.

import os
import logging

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (raw reads)
def get_expected_output_fastqc_raw(species):

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_raw", {}).get("execute", True) == False:
        logging.info(f"Skipping FastQC for raw reads for {species}. Disabled in config.")
        return []

    files = get_r1_read_files_for_species(species)

    all_inputs = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','')
        all_inputs.append(os.path.join(species, "results", "reads", "reads_raw", "fastqc", f"{filename}_raw_fastqc.html"))
        # Add R2 if exists
        if os.path.exists(raw_file.replace("_R1_", "_R2_")):
            all_inputs.append(os.path.join(species, "results", "reads", "reads_raw", "fastqc", f"{filename.replace("_R1_", "_R2_")}_raw_fastqc.html"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (adapter trimmed reads)
def get_expected_output_fastqc_trimmed(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_trimmed", {}).get("execute", True) == False:
        logging.info(f"Skipping FastQC for trimmed reads for {species}. Disabled in config.")
        return []

    all_inputs = []
    for sample in get_sample_ids_for_species(species):
        all_inputs.append(os.path.join(species, "results", "reads", "reads_trimmed", "fastqc",  f"{sample}_trimmed_fastqc.html"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (adapter removed reads)
def get_expected_output_fastqc_quality_filtered(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_quality_filtered", {}).get("execute", True) == False:
        logging.info(f"Skipping FastQC for quality filtered reads for {species}. Disabled in config.")
        return []

    all_inputs = []
    for sample in get_sample_ids_for_species(species):
        all_inputs.append(os.path.join(species, "results", "reads", "reads_quality_filtered", "fastqc", f"{sample}_quality_filtered_fastqc.html"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for FastQC (merged reads)
def get_expected_output_fastqc_merged(species):
    all_inputs = []
    for individual in get_individuals_for_species(species):
        all_inputs.append(os.path.join(species, "results", "reads", "reads_merged", "fastqc", f"{individual}_merged_fastqc.html"))
    return all_inputs

#-----------------------------------------------------------------------------------------------
# Get expected output file paths for contamination analysis (ECMSD)
def get_expected_output_contamination_ecmsd(species):  

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("contamination_analysis", {}).get("tools", {}).get("ecmsd", {}).get("execute", True) == False:
        logging.info(f"Skipping contamination analysis with Centrifuge for {species}. Disabled in config.")
        return []

    expected_outputs = []

    for individual in get_individuals_for_species(species):
        for sample in get_samples_for_species_individual(species, individual):
            expected_outputs.append(os.path.join(species, "results", "contamination_analysis", "ecmsd", individual, sample, "mapping", f"{sample}_Mito_summary.txt"))
    
    return expected_outputs

#-----------------------------------------------------------------------------------------------
# Get expected output file paths for contamination analysis (Centrifuge)
def get_expected_output_contamination_centrifuge(species):  

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("contamination_analysis", {}).get("tools", {}).get("centrifuge", {}).get("execute", True) == False:
        logging.info(f"Skipping contamination analysis with Centrifuge for {species}. Disabled in config.")
        return []

    expected_outputs = []

    for individual in get_individuals_for_species(species):
        for sample in get_samples_for_species_individual(species, individual):
            expected_outputs.append(os.path.join(species, "results", "contamination_analysis", "centrifuge", individual, sample, f"{sample}_centrifuge_report.tsv"))
            expected_outputs.append(os.path.join(species, "results", "contamination_analysis", "centrifuge", individual, sample, f"{sample}_taxon_counts.tsv"))
    
    return expected_outputs

#-----------------------------------------------------------------------------------------------
# Get expected output file paths for contamination analysis (all tools)
def get_expected_output_contamination(species):  

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("contamination_analysis", {}).get("execute", True) == False:
        logging.info(f"Skipping contamination analysis for {species}. Disabled in config.")
        return []

    expected_outputs = []

    # Add ECMSD contamination analysis outputs
    expected_outputs += get_expected_output_contamination_ecmsd(species)

    # Add Centrifuge contamination analysis outputs
    expected_outputs += get_expected_output_contamination_centrifuge(species)

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for MultiQC reports
def get_expected_output_multiqc(species):
    expected_outputs = []

    # Add MultiQC reports for different read processing stages
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_raw", {}).get("execute", True) == True:
        expected_outputs.append(os.path.join(species, "results", "reads", f"{species}_multiqc_raw.html"))
    else:
        logging.info(f"Skipping MultiQC report for raw reads for {species}. Disabled in config.")

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_trimmed", {}).get("execute", True) == True:
        expected_outputs.append(os.path.join(species, "results", "reads", f"{species}_multiqc_trimmed.html"))
    else:
        logging.info(f"Skipping MultiQC report for trimmed reads for {species}. Disabled in config.")

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_quality_filtered", {}).get("execute", True) == True:
        expected_outputs.append(os.path.join(species, "results", "reads", f"{species}_multiqc_quality_filtered.html"))
    else:
        logging.info(f"Skipping MultiQC report for quality filtered reads for {species}. Disabled in config.")
    
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("quality_checking_merged", {}).get("execute", True) == True:
        expected_outputs.append(os.path.join(species, "results", "reads", f"{species}_multiqc_merged.html"))
    else:
        logging.info(f"Skipping MultiQC report for merged reads for {species}. Disabled in config.")

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get expected output file paths for read count plots
def get_expected_output_reads_plots(species):
    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("statistical_analysis", {}).get("execute", True) == False:
        logging.info(f"Skipping read plots generation for {species}. Disabled in config.")
        return []
    
    expected_outputs = []

    expected_outputs.append(os.path.join(species, "results", "reads", "plots", f"{species}_read_counts.png"))
    expected_outputs.append(os.path.join(species, "results", "reads", "plots", f"{species}_read_counts_comparison_by_individual.png"))

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get all expected output file paths for raw read processing
def get_expexted_output_raw_read_processing(species):

    if config.get("pipeline", {}).get("raw_reads_processing", {}).get("execute", True) == False:
        logging.info(f"Skipping raw read processing for {species}. Disabled in config.")
        return []
    
    # Add FastQC outputs for raw reads
    expected_outputs = []

    # Add MultiQC reports for different read processing stages
    expected_outputs += get_expected_output_multiqc(species)

    # Add ECMSD contamination analysis outputs
    expected_outputs += get_expected_output_contamination(species)

    # Add summary plots for read count comparisons
    expected_outputs += get_expected_output_reads_plots(species)

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Get all expected output file paths for reference processing
def get_expected_output_reference_processing(species):

    if config.get("pipeline", {}).get("reference_processing", {}).get("execute", True) == False:
        logging.info(f"Skipping reference processing for {species}. Disabled in config.")
        return []

    expected_outputs = []

    try:
    # Get all reference for the species
        references_list = get_references_ids_for_species(species)

    except Exception as e: 
        # Print error if reference files are missing or inaccessible
        logging.error(e)
        return []
    
     # Get all individuals for the species
    individuals = get_individuals_for_species(species)

    for reference in references_list:

        if config.get("pipeline", {}).get("reference_processing", {}).get("endogenous_reads_analysis", {}).get("execute", True) == True:
            # Add endogenous and coverage plots for each reference and individual
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "endogenous_reads", f"{species}_{reference}_endogenous_reads_bar_chart.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "endogenous_reads", f"{species}_{reference}_raw_and_endogenous_reads_bar_chart.png"))
        else:
            logging.info(f"Skipping endogenous reads analysis for species {species} and reference {reference}. Disabled in config.")
        
        if config.get("pipeline", {}).get("reference_processing", {}).get("coverage_analysis", {}).get("execute", True) == True:
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_depth_coverage_violin.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_depth_coverage_bar.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_coverage_breadth_bar.png"))
            expected_outputs.append(os.path.join(species, "results" ,reference, "plots", "coverage", f"{species}_{reference}_individual_coverage_breadth_violin.png"))

        else:
            logging.info(f"Skipping coverage analysis for species {species} and reference {reference}. Disabled in config.")

        #expected_outputs.append(os.path.join(species, "results" ,reference_id, "multiqc.html"))


        for individual in individuals:
            
            expected_outputs.append(os.path.join(species, "processed" ,reference, "mapped", f"{individual}_{reference}_final.bam"))
            expected_outputs.append(os.path.join(species, "processed" ,reference, "mapped", f"{individual}_{reference}_final.bam.bai"))

            #"{species}/results/{reference}/analytics/{individual}_multiqc.html",
            expected_outputs.append(os.path.join(species, "results" ,reference, "analytics", f"{individual}_{reference}_multiqc.html"))

            if config.get("pipeline", {}).get("reference_processing", {}).get("damage_analysis", {}).get("execute", True) == True:
                #{species}/results/{reference}/analytics/{individual}/mapdamage/misincorporation.txt
                expected_outputs.append(os.path.join(species, "results" ,reference, "analytics", individual, "mapdamage", "misincorporation.txt"))
                #expected_outputs.append(directory(os.path.join(species, "results" ,reference_id, "damage", "damageprofile", ind)))
            else:
                logging.info(f"Skipping damage analysis for species {species} and individual {individual} to reference {reference}. Disabled in config.")

    return expected_outputs

# -----------------------------------------------------------------------------------------------
# Skip existing files based on configuration
def skip_existing_files(expected_outputs):

    # Filter out files that already exist if the skip_existing_files option is enabled
    if config.get("pipeline", {}).get("global", {}).get("skip_existing_files", True) == False:
        return expected_outputs
    
    logging.info("Checking for existing files to skip...")

    expected_outputs_existing = []
    expected_outputs_not_existing = []

    # Check for existing files and filter them out
    for output in expected_outputs:
        if os.path.exists(output):
            expected_outputs_existing.append(output)
        else:
            expected_outputs_not_existing.append(output)

    # Log existing files that will be skipped
    if len(expected_outputs_existing) > 0:
        logging.info("The following files already exist and will be skipped:")
        for existing_file in expected_outputs_existing:
            logging.info("\t" + "- Skipping: " + existing_file)
       
    return expected_outputs_not_existing

# -----------------------------------------------------------------------------------------------
# Generate all input file paths required for the 'all' rule in Snakemake
# This function collects all expected output files for each species and sample, ensuring that
# downstream rules have the correct input targets for completion. It is typically used to define
# the 'all' rule in the Snakefile, which triggers the entire workflow.
def get_expected_outputs_from_pipeline(wildcards):
    # Initialize the list to hold all required input file paths
    expected_output = []

    # Loop over each species defined in the config (must be available in the global scope)
    for species in config["species"]:
        expected_output += get_expexted_output_raw_read_processing(species)
        expected_output += get_expected_output_reference_processing(species)

        expected_output.append(os.path.join(species, "results", "summary", f"{species}_multiqc.overall.html"))

        #{species}/results/analytics/{individual}_multiqc.html
        for individual in get_individuals_for_species(species):
            expected_output.append(os.path.join(species, "results", "summary", f"{individual}_multiqc.html"))

    # Optionally skip files that already exist to avoid redundant processing
    expected_output = skip_existing_files(expected_output)

    # Log all determined inputs for debugging and traceability
    logging.info("Determined input for rule 'all':")
    for input in expected_output:
        logging.info("\t" + "- Requesting: " + input)

    # Return the complete list of input file paths
    return expected_output
