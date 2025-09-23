# =================================================================================================
#     Input Manager Utility Functions for aDNA Pipeline
# =================================================================================================
# This script provides helper functions for managing input files and sample metadata in the pipeline.
# Functions include file discovery, reference genome handling, and sample identification.

import os
import glob
import logging

# -----------------------------------------------------------------------------------------------
# Get all files in a folder matching a specific pattern (e.g., *.fastq.gz)
def get_files_in_folder_matching_pattern(folder: str, pattern: str) -> list:
    # Check if the folder exists
    if not os.path.exists(folder):
        raise Exception(f"Invalid folder: {folder}")
    # Read all files matching the pattern into a list
    files = glob.glob(os.path.join(folder, pattern))
    return files

# -----------------------------------------------------------------------------------------------
# Get reference genome files for a given species (supports .fna, .fasta, .fa)
def get_reference_genome_file_list_for_species(species: str) -> list[tuple[str, str]]:
    # Construct reference genome folder path
    ref_genome_folder = os.path.join(species, "raw", "ref_genome")
    # Collect all supported reference genome files
    reference_genome_files = get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fna")
    reference_genome_files += get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fasta")
    reference_genome_files += get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fa")
    # Raise error if no reference genome files found
    if len(reference_genome_files) == 0:
        raise Exception(f"No reference genome found for species {species}.")
    # Return as list of tuples: (filename without extension, full path)
    reference_genome_files_with_filename = [(os.path.splitext(os.path.basename(f))[0], f) for f in reference_genome_files]
    return reference_genome_files_with_filename

# -----------------------------------------------------------------------------------------------
# Get individual sample IDs for a given species based on raw read filenames
def get_individuals_for_species(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")
    if len(files) == 0:
        raise Exception(f"No raw reads found for species {species}.")
    individuals = set()
    for f in files:
        basename = os.path.basename(f)
        # Extract individual ID (before first underscore)
        individual = basename.split("_")[0]
        individuals.add(individual)
    return sorted(list(individuals))

# -----------------------------------------------------------------------------------------------
# Get only reference genome file paths for a species
def get_ref_genomes_for_species(species):
    refs = get_reference_genome_file_list_for_species(species)
    return [ref[1] for ref in refs]

# -----------------------------------------------------------------------------------------------
# Get sample IDs for a species based on raw read filenames
def get_sample_ids_for_species(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")
    samples = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        samples.append(filename)
    if samples:
        print(f"Found {len(samples)} samples for species {species}.")
        print (f"Samples: {samples}")
    return samples

# -----------------------------------------------------------------------------------------------
# Get input file paths for MultiQC (raw reads)
def get_input_multiqc_raw(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")
    all_inputs = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','')
        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "raw", f"{filename}_fastqc.html"))
        # Add R2 if exists
        if os.path.exists(raw_file.replace("_R1_", "_R2_")):
            all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "raw", f"{filename.replace("_R1_", "_R2_")}_fastqc.html"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get input file paths for reads processing statistics (raw)
def get_input_reads_processing_raw(species):
    all_inputs = []
    for sample in get_sample_ids_for_species(species):
        all_inputs.append(os.path.join(species, "processed", "qualitycontrol", "statistics", f"{sample}_all_steps_count.csv"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get input file paths for trimmed reads statistics
def get_input_reads_processing_trimmed(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")
    all_inputs = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        all_inputs.append(os.path.join(species, "processed", "qualitycontrol", "statistics", f"{filename}_trimmed_reads.count"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get input file paths for quality filtered reads statistics
def get_input_reads_processing_quality_filtered(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")
    all_inputs = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        all_inputs.append(os.path.join(species, "processed", "qualitycontrol", "statistics", f"{filename}_quality_filtered_reads.count"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get input file paths for MultiQC (trimmed reads)
def get_input_multiqc_trimmed(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")
    all_inputs = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "quality_filtered", f"{filename}_quality_filtered_fastqc.html"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get input file paths for MultiQC (quality filtered reads)
def get_input_multiqc_quality_filtered(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")
    all_inputs = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "quality_filtered", f"{filename}_quality_filtered_fastqc.html"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Get input file paths for MultiQC (merged reads)
def get_input_multiqc_merged(species):
    individuals = get_individuals_for_species(species)
    all_inputs = []
    for individual in individuals:
        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "merged", f"{individual}_fastqc.html"))
    return all_inputs

# -----------------------------------------------------------------------------------------------
# Generate all input file paths required for the 'all' rule in Snakemake
# This function collects all expected output files for each species and sample, ensuring that
# downstream rules have the correct input targets for completion. It is typically used to define
# the 'all' rule in the Snakefile, which triggers the entire workflow.
def get_all_inputs(wildcards):
    # Initialize the list to hold all required input file paths
    all_inputs = []
    # Loop over each species defined in the config (must be available in the global scope)
    for species in config["species"]:
        species_folder = species
        # Add the main quality check report for each species
        all_inputs.append(os.path.join(species_folder, "results","qualitycontrol", f"quality_check_report_{species}.html"))
        # For each sample, add the ECMSD mitochondrial summary file
        for sample in get_sample_ids_for_species(species):
            all_inputs.append(os.path.join(species_folder, "results", "qualitycontrol", "ecmsd", sample, "mapping", "Mito_summary.txt"))
        
        # Add summary plots for read count comparisons
        all_inputs.append(os.path.join(species_folder, "results", "plots", f"{species}_read_count_comparison.png"))
        all_inputs.append(os.path.join(species_folder, "results", "plots", f"{species}_read_count_comparison_by_individual.png"))
        
        # Get all individuals for the species
        individuals = get_individuals_for_species(species)

        try:
            # Get all reference genomes for the species
            ref_genome_list = get_reference_genome_file_list_for_species(species)
            for ref_genome_tuple in ref_genome_list:

                ref_genome_id = ref_genome_tuple[0]
                # Add endogenous and coverage plots for each reference genome and individual
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "endogenous_reads", f"{species}_{ref_genome_id}_endogenous_reads_pie_chart.pdf"))
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "endogenous_reads", f"{species}_{ref_genome_id}_endogenous_reads_bar_chart.png"))
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_breadth_coverage_violin.png"))
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_depth_coverage_violin.png"))
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_depth_coverage_bar.png"))
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_breadth_coverage_bins.png"))
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_depth_violin.png"))                
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_coverage_breadth_bar.png"))
                all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_coverage_breadth_violin.png"))
                # Add consensus sequence output for each individual and reference genome
                all_inputs.append(os.path.join(species_folder, "processed" ,ref_genome_id, "consensus", f"{ind}_{ref_genome_id}", f"{ind}_{ref_genome_id}_consensus.fa.gz"))
        except Exception as e: 
            # Print error if reference genome files are missing or inaccessible
            print(e)
            pass

        #for ind in individuals:
            
    # Log all determined inputs for debugging and traceability
    logging.info("Determined input for rule 'all':")
    for input in all_inputs:
        logging.info("\t" + "- " + input)
    # Return the complete list of input file paths
    return all_inputs
