import os
from common_aDNA_scripts import *

def fastqc_for_raw_data(species: str):
    print_species_execution(f"Running fastqc for species {species} raw data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_DUPLICATES_REMOVED,
        species
    ):
        print_skipping(f"FASTQC for raw is not enabled for {species} in the config. Skipping this step.")
        return

    # raw data
    raw_reads_folder = get_folder_path_species_raw_reads(species)
    output_folder = get_folder_path_species_results_qc_fastqc_raw(species)

    raw_reads_files = get_files_in_folder_matching_pattern(raw_reads_folder, f"*{FILE_ENDING_FASTQ_GZ}")

    files_list = get_files_in_folder_matching_pattern(output_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(files_list) > 0:
        print_skipping(f"Fastqc for raw already exists for species {species}.")
        return
    
    if len(raw_reads_files) == 0:
        print_warning(f"No raw reads found for species {species}.")
        return
    
    execute_fastqc(species, raw_reads_files, output_folder)

    print_success(f"fastqc for species {species} raw data complete")

def fastqc_for_quality_filtered_data(species: str):
    print_species_execution(f"Running fastqc for species {species} quality filtered data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_DUPLICATES_REMOVED,
        species
    ):
        print_skipping(f"FASTQC for quality filtered is not enabled for {species} in the config. Skipping this step.")
        return

    quality_filtered_reads_folder = get_folder_path_species_processed_quality_filtered(species)
    output_folder = get_folder_path_species_results_qc_fastqc_quality_filtered(species)

    quality_filtered_reads = get_files_in_folder_matching_pattern(quality_filtered_reads_folder, f"*{FILE_ENDING_QUALITY_FILTERED_FASTQ_GZ}")

    files_list = get_files_in_folder_matching_pattern(output_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(files_list) > 0:
        print_skipping(f"Fastqc for quality filtered already exists for species {species}.")
        return
    
    if len(quality_filtered_reads) == 0:
        print_warning(f"No quality filtered reads found for species {species}.")
        return
    
    execute_fastqc(species, quality_filtered_reads, output_folder)

    print_success(f"fastqc for species {species} quality filtered data complete")

def fastqc_for_adapter_removed_data(species: str):
    print_species_execution(f"Running fastqc for species {species} adapter removed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_DUPLICATES_REMOVED,
        species
    ):
        print_skipping(f"FASTQC for adapter removed is not enabled for {species} in the config. Skipping this step.")
        return

    #adapter removed data
    trimmed_reads_folder = get_folder_path_species_processed_adapter_removed(species)
    output_folder = get_folder_path_species_results_qc_fastqc_adapter_removed(species)

    adapter_removed_reads = get_files_in_folder_matching_pattern(trimmed_reads_folder, f"*{FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ}")

    files_list = get_files_in_folder_matching_pattern(output_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(files_list) > 0:
        print_skipping(f"fastqc for adapter removed already exists for species {species}.")
        return
    
    if len(adapter_removed_reads) == 0:
        print_warning(f"No adapter removed reads found for species {species}.")
        return
    
    execute_fastqc(species, adapter_removed_reads, output_folder)
        
    print_success(f"fastqc for species {species} adapter removed data complete")

def fastqc_for_duplicates_removed_data(species: str):
    print_species_execution(f"Running fastqc for species {species} duplicates removed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_DUPLICATES_REMOVED,
        species
    ):
        print_skipping(f"FASTQC for duplicates removed is not enabled for {species} in the config. Skipping this step.")
        return

    duplicates_removed_reads_folder = get_folder_path_species_processed_duplicates_removed(species)
    output_folder = get_folder_path_species_results_qc_fastqc_duplicates_removed(species)

    duplicates_removed_reads = get_files_in_folder_matching_pattern(duplicates_removed_reads_folder, f"*{FILE_ENDING_DUPLICATES_REMOVED_FASTQ_GZ}")

    files_list = get_files_in_folder_matching_pattern(output_folder, f"*{FILE_ENDING_FASTQC_HTML}")

    if len(files_list) > 0:
        print_skipping(f"Fastqc for duplicates removed already exists for species {species}.")
        return
    
    if len(duplicates_removed_reads) == 0:
        print_warning(f"No duplicate removed reads found for species {species}.")
        return
    
    execute_fastqc(species, duplicates_removed_reads, output_folder)

    print_success(f"fastqc for species {species} duplicates removed data complete")

def execute_fastqc(species: str, reads_file_list: list, output_folder: str):

    if len(reads_file_list) == 0:
        raise Exception(f"No reads provided.")

    print_info(f"Running fastqc for {len(reads_file_list)} files for species {species}")

    threads = config.get_threads()
    

    command = [
        PROGRAM_PATH_FASTQC,
        "-o", output_folder,
        "-t", str(threads),
        *reads_file_list
    ]

    try:
        run_command(
            command, 
            description=f"FastQC for species {species}"
        )
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to run FastQC for species {species}: {e}")


def all_species_fastqc_raw():
    print_step_execution("Running fastqc for all species raw data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_RAW
    ):
        print_skipping(f"FASTQC for raw data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        fastqc_for_raw_data(species)

    print_info("Fastqc for all species raw data completed successfully.")

def all_species_fastqc_adapter_removed():
    print_step_execution("Running fastqc for all species trimmed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_ADAPTER_REMOVED
    ):
        print_skipping(f"FASTQC for adapter removed data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
       fastqc_for_adapter_removed_data(species)

    print_info("Fastqc for all species trimmed data completed successfully.")

def all_species_fastqc_quality_filtered():
    print_step_execution("Running fastqc for all species quality filtered data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_QUALITY_FILTERED
    ):
        print_skipping(f"FASTQC for quality filtered data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        fastqc_for_quality_filtered_data(species)

    print_info("Fastqc for all species quality filtered data completed successfully.")

def all_species_fastqc_duplicates_removed():
    print_step_execution("Running fastqc for all species duplicates removed data")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_CHECK_RAW
    ):
        print_skipping(f"FASTQC for duplicates removed data is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        fastqc_for_duplicates_removed_data(species)

    print_info("Fastqc for all species duplicates removed data completed successfully.")