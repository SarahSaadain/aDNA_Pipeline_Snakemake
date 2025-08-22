import os
from common_aDNA_scripts import *

import raw_reads_processing.common_raw_reads_processing_helpers as common_rrp

def execute_fastp_deduplication(input_file_path:str, output_file_path:str, threads:int = THREADS_DEFAULT):
    print_info(f"Filtering {input_file_path} ...")

    if not os.path.exists(input_file_path):
        raise Exception(f"Read file {input_file_path} does not exist!")
    
    if os.path.exists(output_file_path):
        print_skipping(f"Output file {output_file_path} already exists!")
        return
    
    filepath_reads_failed = output_file_path.replace(FILE_ENDING_DUPLICATES_REMOVED_FASTQ_GZ, "_failed.fastq.gz")
    filepath_json_report = output_file_path.replace(FILE_ENDING_DUPLICATES_REMOVED_FASTQ_GZ, FILE_ENDING_FASTP_JSON_REPORT)
    filepath_html_report = output_file_path.replace(FILE_ENDING_DUPLICATES_REMOVED_FASTQ_GZ, FILE_ENDING_FASTP_HTML_REPORT)

    #https://github.com/OpenGene/fastp/blob/59cc2f67414e74e99d42774e227b192a3d9bb63a/README.md#all-options
    command_fastp = [
        PROGRAM_PATH_FASTP, 
        "--dedup",                              #enable deduplication to drop the duplicated reads/pairs
        "--disable_adapter_trimming",
        "--disable_length_filtering",
        "--disable_quality_filtering",
        "--thread", str(threads),                    # Number of threads
        "--in1", input_file_path,               # Input R1 file
        "--out1", output_file_path,
        "--failed_out", filepath_reads_failed,
        "--json", filepath_json_report,
        "--html", filepath_html_report
    ]

    print_debug(f"Executing command: {' '.join(command_fastp)}")
    
    try:
        run_command(
            command_fastp, 
            description=f"fastp deduplication for {input_file_path}"
        )
        print_success(f"fastp deduplication for {input_file_path} complete")
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to run fastp deduplication for {input_file_path}: {e}")


def fastp_deduplication_for_species(species: str):

    print_species_execution(f"Running fastp deduplication for {species}")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.DEDUPLICATION,
        species
    ):
        print_skipping(f"Deduplication is not enabled for {species} in the config. Skipping this step.")
        return

    reads_folder = get_folder_path_species_processed_quality_filtered(species)

    reads_files_list = get_files_in_folder_matching_pattern(reads_folder, f"*{FILE_ENDING_QUALITY_FILTERED_FASTQ_GZ}")

    number_of_reads_files = len(reads_files_list)

    if number_of_reads_files == 0:
        print_warning(f"No quality filtered reads found for species {species}.")
        return
    
    print_debug(f"Found {len(reads_files_list)} quality filtered reads for species {species}.")
    print_debug(f"Quality filtered reads: {reads_files_list}")

    count_current = 0
    
    for read_file_path in reads_files_list:

        count_current += 1

        print_info(f"[{count_current}/{number_of_reads_files}] Processing read file: {get_filename_from_path(read_file_path)}")
        
        individual = common_rrp.get_individual_from_file(read_file_path)

        is_ref_genome_read_file_exists = common_rrp.is_species_individual_reads_file_exists(
            species, 
            individual
        )

        if is_ref_genome_read_file_exists:
            print_skipping(f"Individual {individual} already prepared for reference genome processing!")
            continue

        output_file_path = common_rrp.get_deduplication_path_for_quality_filtered_reads(species, read_file_path)
        execute_fastp_deduplication(read_file_path, output_file_path)

    print_info(f"fastp deduplication for {species} complete")

def all_species_fastp_deduplication():
    print_step_execution("Running fastp deduplication for all species")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.DEDUPLICATION
    ):
        print_skipping("Deduplication is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        fastp_deduplication_for_species(species)

    print_info("fastp deduplication for all species complete")
