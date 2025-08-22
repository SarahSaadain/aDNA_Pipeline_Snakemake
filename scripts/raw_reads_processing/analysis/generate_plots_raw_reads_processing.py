import os
from common_aDNA_scripts import *

# Function to plot the number of reads before and after processing for a single species.
def plot_reads_processing_result(species: str):
    print_info(f"Plotting reads processing result for species {species}")

    # Build the path to the TSV file that contains reads processing results.
    input_file_path = os.path.join(
        get_folder_path_species_results_qc_reads_processing(species),
        f"{species}{FILE_ENDING_READS_PROCESSING_RESULT_TSV}"
    )

    # Check if the input file exists before trying to use it.
    if not os.path.exists(input_file_path):
        print_warning(f"Input file not found: {input_file_path}")
        return  # Skip plotting if file doesn't exist

    # Define where the output plot should be saved.
    output_folder_path = get_folder_path_species_results_plots_reads_processing(species)

    # Retrieve the path to the R script that plots reads before/after processing.
    r_script = get_r_script(R_SCRIPT_PLOT_READS_BEFORE_AFTER_PROCESSING, FOLDER_RAW_READS_PROCESSING)

    # Run the R script with species name, input TSV path, and output folder.
    call_r_script(r_script, species, input_file_path, output_folder_path)

# Function to generate sequence length distribution plots for a given species.
def plot_read_length_distribution(species: str):
    print_info(f"Plotting read length distribution for species {species}")

    # Get all TSV files in the read length distribution folder for this species.
    analysis_files = get_files_in_folder_matching_pattern(
        get_folder_path_species_results_qc_read_length_distribution(species),
        f"*{FILE_ENDING_READ_LENGTH_DISTRIBUTION_TSV}"
    )

    # If no matching files are found, skip this step.
    if len(analysis_files) == 0:
        print_warning(f"No read length distribution files found for species {species}.")
        return

    # Retrieve the R script that generates the length distribution plots.
    r_script = get_r_script(R_SCRIPT_PLOT_SEQUENCE_LENGTH_DISTRIBUTION, FOLDER_RAW_READS_PROCESSING)

    # Loop through all files and plot each one.
    for analysis_file in analysis_files:
        # Define output folder for the plot.
        output_folder_path = get_folder_path_species_results_plots_read_length_distribution(species)

        print_info(f"Plotting read length distribution for file {analysis_file} to {output_folder_path}")

        # Run the R script with species, file, and output folder.
        call_r_script(r_script, species, analysis_file, output_folder_path)

    print_info(f"Finished plotting read length distribution for species {species}")

def  plot_contamination(species):
    print_info(f"Plotting contamination analysis for species {species}")

    # Get the path to the contamination analysis file.
    kraken_contamination_by_invdividual = os.path.join(
        get_folder_path_species_results_qc_kraken(species),
        f"{species}{FILE_ENDING_KRAKEN_BY_INDIVIDUAL_COMBINED_ANALYSIS_CSV}"
    )

    # Check if the input file exists before trying to use it.
    if not os.path.exists(kraken_contamination_by_invdividual):
        print_warning(f"Input file not found: {kraken_contamination_by_invdividual}")
        return  # Skip plotting if file doesn't exist

    # check if file is empty
    if os.path.getsize(kraken_contamination_by_invdividual) == 0:
        print_warning(f"Input file is empty: {kraken_contamination_by_invdividual}")
        return

    # Define where the output plot should be saved.
    output_folder_path = get_folder_path_species_results_plots_contamination(species)

    # Retrieve the R script for plotting contamination analysis.
    r_script = get_r_script(R_SCRIPT_PLOT_CONTAMINATION_KRAKEN, FOLDER_RAW_READS_PROCESSING)

    # Run the R script with species name, input TSV path, and output folder.
    call_r_script(r_script, species, kraken_contamination_by_invdividual, output_folder_path)

# Master function to run all plotting routines for a single species.
def species_generate_plots(species: str):
    print_species_execution(f"Generating plots for species {species}")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.GENERATE_RAW_READS_PLOTS,
        species
    ):
        print_skipping(f"Generate raw reads plots is not enabled for {species} in the config. Skipping this step.")
        return

    # Plot reads before/after processing.
    plot_reads_processing_result(species)

    # Plot sequence length distribution.
    plot_read_length_distribution(species)

    # Plot contamination analysis.
    plot_contamination(species)

    print_info(f"Finished generating plots for species {species}")

# Run all plotting functions for every species defined in the configuration.
def all_species_generate_plots():
    print_step_execution("Generating plots for all species")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.GENERATE_RAW_READS_PLOTS
    ):
        print_skipping(f"Generate raw reads plots is not enabled in the config. Skipping this step.")
        return

    # Loop through each species and generate plots.
    for species in FOLDER_SPECIES:
        species_generate_plots(species)

    print_info("Finished generating plots for all species")