from common_aDNA_scripts import *

def multiqc_for_raw_data(species: str):
    print_species_execution(f"Running MultiQC for species {species} raw data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_RAW,
        species
    ):
        print_skipping(f"MULTIQC for raw data for {species} is not enabled in the config. Skipping this step.")
        return

    # raw data
    raw_fastqc_folder = get_folder_path_species_results_qc_fastqc_raw(species)

    file_list = get_files_in_folder_matching_pattern(raw_fastqc_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(file_list) == 0:
        print_warning(f"No fastqc data found for species {species}.")
        return

    output_folder = get_folder_path_species_results_qc_multiqc_raw(species)

    files_list = get_files_in_folder_matching_pattern(output_folder, "*.html")

    if len(files_list) > 0:
        print_skipping(f"Multiqc for raw already exists for species {species}.")
        return

    run_multiqc(species, raw_fastqc_folder, output_folder)

    print_success(f"MultiQC for species {species} raw data complete")

def multiqc_for_quality_filtered_data(species: str):
    print_species_execution(f"Running MultiQC for species {species} quality filtered data")  

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_QUALITY_FILTERED,
        species
    ):
        print_skipping(f"MULTIQC for quality filtered data for {species} is not enabled in the config. Skipping this step.")
        return

    quality_filtered_fastqc_folder = get_folder_path_species_results_qc_fastqc_quality_filtered(species)

    file_list = get_files_in_folder_matching_pattern(quality_filtered_fastqc_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(file_list) == 0:
        print_warning(f"No fastqc data found for species {species}.")
        return

    output_folder = get_folder_path_species_results_qc_multiqc_quality_filtered(species)

    files_list_multiqc_check = get_files_in_folder_matching_pattern(output_folder, "*.html")

    if len(files_list_multiqc_check) > 0:
        print_skipping(f"Multiqc for quality filtered already exists for species {species}.")
        return

    run_multiqc(species, quality_filtered_fastqc_folder, output_folder)

    print_success(f"MultiQC for species {species} quality filtered data complete")

def multiqc_for_duplicates_removed_data(species: str):
    print_species_execution(f"Running MultiQC for species {species} duplicates removed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_DUPLICATES_REMOVED,
        species
    ):
        print_skipping(f"MULTIQC for duplicates removed data for {species} is not enabled in the config. Skipping this step.")
        return

    duplicates_removed_fastqc_folder = get_folder_path_species_results_qc_fastqc_duplicates_removed(species)

    file_list = get_files_in_folder_matching_pattern(duplicates_removed_fastqc_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(file_list) == 0:
        print_warning(f"No fastqc data found for species {species}.")
        return

    output_folder = get_folder_path_species_results_qc_multiqc_duplicates_removed(species)

    files_list = get_files_in_folder_matching_pattern(output_folder, "*.html")

    if len(files_list) > 0:
        print_skipping(f"Multiqc for duplicates already exists for species {species}.")
        return

    run_multiqc(species, duplicates_removed_fastqc_folder, output_folder)

    print_success(f"MultiQC for species {species} duplicates removed data complete")

def multiqc_for_adapter_removed_data(species: str):
    print_species_execution(f"Running MultiQC for species {species} trimmed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_ADAPTER_REMOVED,
        species
    ):
        print_skipping(f"MULTIQC for adapter removed data for {species} is not enabled in the config. Skipping this step.")
        return

    #adapter removed data
    trimmed_fastqc_folder = get_folder_path_species_results_qc_fastqc_adapter_removed(species)

    file_list = get_files_in_folder_matching_pattern(trimmed_fastqc_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(file_list) == 0:
        print_warning(f"No fastqc data found for species {species}.")
        return

    output_folder = get_folder_path_species_results_qc_multiqc_adapter_removed(species)

    files_list = get_files_in_folder_matching_pattern(output_folder, "*.html")

    if len(files_list) > 0:
        print_skipping(f"Multiqc for adapter removed already exists for species {species}.")
        return

    run_multiqc(species, trimmed_fastqc_folder, output_folder)

    print_success(f"Multiqc for species {species} trimmed data complete")

def run_multiqc(species: str, fastqc_results_folder: str, output_folder: str):

    command = [PROGRAM_PATH_MULTIQC, fastqc_results_folder, "-o", output_folder]

    try:
        run_command(
            command, 
            description=f"MultiQC for species {species}"
        )
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to run MultiQC for species {species}: {e}")

def all_species_multiqc_raw():
    print_step_execution("Running MultiQC for all species raw data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_RAW
    ):
        print_skipping(f"MULTIQC for raw data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        multiqc_for_raw_data(species)

    print_success("Multiqc for all species raw data completed successfully.")

def all_species_multiqc_adapter_removed():
    print_step_execution("Running MultiQC for all species trimmed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_ADAPTER_REMOVED
    ):
        print_skipping(f"MULTIQC for adapter removed data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
       multiqc_for_adapter_removed_data(species)

    print_success("Multiqc for all species trimmed data completed successfully.")

def all_species_multiqc_quality_filtered():
    print_step_execution("Running MultiQC for all species quality filtered data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_QUALITY_FILTERED
    ):
        print_skipping(f"MULTIQC for quality filtered data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        multiqc_for_quality_filtered_data(species)

    print_success("Multiqc for all species quality filtered data completed successfully.")

def all_species_multiqc_duplicates_removed():
    print_step_execution("Running MultiQC for all species duplicates removed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_DUPLICATES_REMOVED
    ):
        print_skipping(f"MULTIQC for duplicates removed data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        multiqc_for_duplicates_removed_data(species)

    print_success("Multiqc for all species duplicates removed data completed successfully.")