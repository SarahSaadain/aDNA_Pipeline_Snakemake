import os
import common_aDNA_scripts as common_adna

# This function runs an R script that generates comparison plots
# showing the number of reads before and after processing for each species.
def plot_comparison_reads_processing_results():
    # Build the path to the folder where R scripts for species comparisons are stored.
    scripts_base_folder = os.path.join(common_adna.FOLDER_ADDITIONAL_ANALYSIS, common_adna.FOLDER_SPECIES_COMPARISON)

    common_adna.print_info(f"Plotting reads before and after comparison")

    # Define the root folder of the aDNA project and the output folder for plots.
    adna_project_folder_path = common_adna.PATH_ADNA_PROJECT
    output_folder_path_for_plot = common_adna.get_folder_path_results_plots()

    # Construct the full path to the specific R script to be used for this plot.
    r_script = common_adna.get_r_script(common_adna.R_SCRIPT_PLOT_COMPARE_SPECIES_READS_BEFORE_AFTER_PROCESSING, scripts_base_folder)

    # Call the R script with the required arguments (project root, config, output folder).
    common_adna.call_r_script(r_script, adna_project_folder_path, common_adna.config_file_path, output_folder_path_for_plot)

    common_adna.print_info(f"Finished reads before and after comparison")

# This function runs an R script to generate plots comparing sequencing depth and breadth
# across multiple species as defined in the configuration.
def plot_comparison_depth_breadth_analysis():
    # Build the base folder path for R scripts related to comparisons.
    scripts_base_folder = os.path.join(common_adna.FOLDER_ADDITIONAL_ANALYSIS, common_adna.FOLDER_SPECIES_COMPARISON)

    common_adna.print_info(f"Plotting depth and breadth comparison.")

    # Define root and output folders.
    adna_project_folder_path = common_adna.PATH_ADNA_PROJECT
    output_folder_path_for_plot = common_adna.get_folder_path_results_plots()

    # Get the R script that performs depth/breadth comparison.
    r_script = common_adna.get_r_script(common_adna.R_SCRIPT_PLOT_COMPARE_SPECIES_DEPTH_BREADTH, scripts_base_folder)

    # Execute the R script.
    common_adna.call_r_script(r_script, adna_project_folder_path, common_adna.config_file_path, output_folder_path_for_plot)

    common_adna.print_info(f"Finished plotting depth and breadth comparison.")

# This function runs an R script that plots comparisons of endogenous read counts
# (i.e., reads mapped to reference genome) across species.
def plot_comparison_endogenous_reads():
    # Set path to the scripts directory for species comparisons.
    scripts_base_folder = os.path.join(common_adna.FOLDER_ADDITIONAL_ANALYSIS, common_adna.FOLDER_SPECIES_COMPARISON)

    common_adna.print_info(f"Plotting endogenous reads comparison.")

    # Set root and output directories.
    adna_project_folder_path = common_adna.PATH_ADNA_PROJECT
    output_folder_path_for_plot = common_adna.get_folder_path_results_plots()

    # Retrieve the R script that compares endogenous reads.
    r_script = common_adna.get_r_script(common_adna.R_SCRIPT_PLOT_COMPARE_SPECIES_ENDOGENOUS_READS, scripts_base_folder)

    # Call the R script with required parameters.
    common_adna.call_r_script(r_script, adna_project_folder_path, common_adna.config_file_path, output_folder_path_for_plot)

    common_adna.print_info(f"Finished plotting endogenous reads comparison.")

# Master function that runs all three comparison plot scripts.
def species_generate_comparison_plots():

    if not common_adna.config.is_process_step_enabled(
        common_adna.PipelineStages.POST_PROCESSING, 
        common_adna.PostProcessingSteps.MTDNA_ANALYSIS.value
    ):
        common_adna.print_skipping("mapDamage analysis is not enabled in the config. Skipping this step.")
        return

    common_adna.print_step_execution("Generating comparison plots for species based on config")

    # Run plotting functions for each analysis type.
    plot_comparison_reads_processing_results()
    plot_comparison_depth_breadth_analysis()
    plot_comparison_endogenous_reads()

    common_adna.print_info(f"Finished generating comparison plots")