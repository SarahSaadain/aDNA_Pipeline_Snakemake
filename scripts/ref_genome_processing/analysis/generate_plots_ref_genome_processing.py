import os
from common_aDNA_scripts import *
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

def plot_depth_analysis(species: str, reference_genome_id: str):
    print_info(f"Plotting depth analysis for species {species} and reference genome {reference_genome_id}")

    analysis_folder = get_folder_path_species_results_refgenome_coverage(species, reference_genome_id)

    print_debug(f"Analysis folder: {analysis_folder}")
    print_debug(f"looking for files with pattern *{FILE_ENDING_COMBINED_COVERAGE_ANALYSIS_CSV}")

    # gives a list with the path and file names
    analysis_files = get_files_in_folder_matching_pattern(analysis_folder, f"*{FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV}")

    # here have multiple files, one for each sample, hence the list
    if len(analysis_files) == 0:
        print_warning(f"No depth analysis files found for species {species} and reference genome {reference_genome_id}.")
        return
    
    print_debug(f"Found {len(analysis_files)} depth analysis files for species {species} and reference genome {reference_genome_id}.")
    print_debug(f"Depth analysis files: {analysis_files}")

    r_script = get_r_script(R_SCRIPT_PLOT_DEPTH, FOLDER_REF_GENOME_PROCESSING)

    print_info(f"R script for plotting depth analysis per sample")

    for analysis_file in analysis_files:

        sample = get_filename_from_path(analysis_file).replace(FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV, "")

        output_folder_path = get_folder_path_species_results_refgenome_plots_depth_sample(species, reference_genome_id, sample)

        print_info(f"Plotting depth analysis for file {analysis_file} to {output_folder_path}")

        call_r_script(r_script, species, analysis_file, output_folder_path)

    
    # plot individuals together
    r_script_individuals = get_r_script(R_SCRIPT_PLOT_DEPTH_COMPARE_INDIVIDUALS, FOLDER_REF_GENOME_PROCESSING)
    output_folder_path_individuals = get_folder_path_species_results_refgenome_plots_depth(species, reference_genome_id)
    combined_results_file = os.path.join(analysis_folder, f"{species}{FILE_ENDING_COMBINED_COVERAGE_ANALYSIS_DETAILED_CSV}")

    print_info(f"Plotting depth analysis for all individuals in {combined_results_file} together to {output_folder_path_individuals}")
    call_r_script(r_script_individuals, species, combined_results_file, output_folder_path_individuals)

    print_info(f"Finished plotting depth analysis for species {species} and reference genome {reference_genome_id}")

def plot_breadth_analysis(species: str, reference_genome_id: str):
    print_info(f"Plotting breadth analysis for species {species}")

    analysis_folder = get_folder_path_species_results_refgenome_coverage(species, reference_genome_id)

    print_debug(f"Analysis folder: {analysis_folder}")
    print_debug(f"looking for files with pattern *{FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV}")

    analysis_files = get_files_in_folder_matching_pattern(analysis_folder, f"*{FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV}")

    if len(analysis_files) == 0:
        print_warning(f"No breadth analysis files found for species {species}.")
        return
    
    print_debug(f"Found {len(analysis_files)} depth analysis files for species {species} and reference genome {reference_genome_id}.")
    print_debug(f"Depth analysis files: {analysis_files}")

    r_script = get_r_script(R_SCRIPT_PLOT_BREADTH, FOLDER_REF_GENOME_PROCESSING)

    print_info(f"R script for plotting breadth analysis per sample")

    for analysis_file in analysis_files:

        sample = get_filename_from_path(analysis_file).replace(FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV, "")

        output_folder_path = get_folder_path_species_results_refgenome_plots_breadth_sample(species, reference_genome_id, sample)

        print_info(f"Plotting breadth analysis for file {analysis_file} to {output_folder_path}")

        call_r_script(r_script, species, analysis_file, output_folder_path)

    # plot individuals together
    r_script_individuals = get_r_script(R_SCRIPT_PLOT_BREADTH_COMPARE_INDIVIDUALS, FOLDER_REF_GENOME_PROCESSING)
    output_folder_path_individuals = get_folder_path_species_results_refgenome_plots_breadth(species, reference_genome_id)

    combined_results_file = os.path.join(analysis_folder, f"{species}{FILE_ENDING_COMBINED_COVERAGE_ANALYSIS_DETAILED_CSV}")

    print_info(f"Plotting breadth analysis for all individuals in {combined_results_file} together to {output_folder_path_individuals}")
    call_r_script(r_script_individuals, species, combined_results_file, output_folder_path_individuals)

    print_info(f"Finished plotting breadth analysis for species {species} and reference genome {reference_genome_id}")

def plot_endogenous_reads(species: str, reference_genome_id: str):
    print_info(f"Plotting endogenous reads for species {species} and reference genome {reference_genome_id}")

    endogenous_reads_analysis_folder = get_folder_path_species_results_refgenome_endogenous_reads(species, reference_genome_id)

    print_debug(f"Analysis folder: {endogenous_reads_analysis_folder}")
    print_debug(f"looking for files with pattern *{FILE_ENDING_ENDOGENOUS_READS_CSV}")

    analysis_files = get_files_in_folder_matching_pattern(endogenous_reads_analysis_folder, f"*{FILE_ENDING_ENDOGENOUS_READS_CSV}")

    if len(analysis_files) == 0:
        print_warning(f"No endogenous reads files found for species {species} and reference genome {reference_genome_id}.")
        return
    
    print_debug(f"Found {len(analysis_files)} depth analysis files for species {species} and reference genome {reference_genome_id}.")
    print_debug(f"Depth analysis files: {analysis_files}")

    r_script = get_r_script(R_SCRIPT_PLOT_ENDOGENOUS_READS, FOLDER_REF_GENOME_PROCESSING)

    for analysis_file in analysis_files:

        output_folder_path = get_folder_path_species_results_refgenome_plots_endogenous_reads(species, reference_genome_id)

        print_info(f"Plotting endogenous reads for file {analysis_file} to {output_folder_path}")

        call_r_script(r_script, species, analysis_file, output_folder_path)

    print_info(f"Finished plotting endogenous reads for species {species} and reference genome {reference_genome_id}")

def plot_snp_data(species: str, reference_genome_id: str):
    print_info(f"Plotting snp data for species {species} and reference genome {reference_genome_id}")

    analysis_folder = get_folder_path_species_processed_refgenome_consensus_sequences(species, reference_genome_id)

    print_debug(f"Analysis folder: {analysis_folder}")
    print_debug(f"looking for files with pattern *{FILE_ENDING_ENDOGENOUS_READS_CSV}")

    analysis_files = get_files_in_folder_matching_pattern(analysis_folder, f"*.mafs.gz")

    if len(analysis_files) == 0:
        print_warning(f"No mafs files found for species {species} and reference genome {reference_genome_id}.")
        return
    
    print_debug(f"Found {len(analysis_files)} mafs files for species {species} and reference genome {reference_genome_id}.")
    print_debug(f"Mafs files: {analysis_files}")

    r_script = get_r_script(R_SCRIPT_PLOT_SNP_BY_INDIVIDUAL, FOLDER_REF_GENOME_PROCESSING)

    output_folder_path = get_folder_path_species_results_refgenome_plots_snp(species, reference_genome_id)

    for analysis_file in analysis_files:

        print_info(f"Plotting snp data for file {analysis_file} to {output_folder_path}")

        call_r_script(r_script, species, analysis_file, output_folder_path)

    r_script = get_r_script(R_SCRIPT_PLOT_SNP_COMPARE_INDIVIDUALS, FOLDER_REF_GENOME_PROCESSING)

    call_r_script(r_script, species, analysis_folder, output_folder_path)

    print_info(f"Finished plotting snp data for species {species} and reference genome {reference_genome_id}")

def species_generate_plots(species: str):
    print_species_execution(f"Generating reference genome plots for species {species}")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.GENERATE_REF_GENOME_PLOTS,
        species
    ):
        print_skipping(f"Generating reference genome plots is not enabled in the config for species {species}. Skipping this step.")
        return

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_tuple in ref_genome_list:

        # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
        ref_genome_id = ref_genome_tuple[0]
        #ref_genome_path = ref_genome_tuple[1]

        print_info(f"Generating plots for reference genome {ref_genome_id}")

        plot_depth_analysis(species, ref_genome_id)
        plot_breadth_analysis(species, ref_genome_id)
        plot_endogenous_reads(species, ref_genome_id)
        plot_snp_data(species, ref_genome_id)

    print_info(f"Finished generating reference genome plots for species {species}")

def all_species_generate_plots():

    print_step_execution("Generating reference genome plots for all species")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.GENERATE_REF_GENOME_PLOTS
    ):
        print_skipping("Generating reference genome plots is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES:
        species_generate_plots(species)

    print_success("Finished generating reference genome plots for all species")
