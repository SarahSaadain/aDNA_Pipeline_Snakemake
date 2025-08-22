import os
from common_aDNA_scripts import *

import raw_reads_processing.common_raw_reads_processing_helpers as common_rrp

def execute_fastp_quality_filter(input_file_path:str, output_file_path:str, threads:int = THREADS_DEFAULT):
    print_info(f"Filtering {input_file_path} ...")

    if not os.path.exists(input_file_path):
        raise Exception(f"Read file {input_file_path} does not exist!")
    
    if os.path.exists(output_file_path):
        print_skipping(f"Output file {output_file_path} already exists!")
        return
    
    filepath_failed_reads = output_file_path.replace(FILE_ENDING_QUALITY_FILTERED_FASTQ_GZ, "_failed.fastq.gz")
    filepath_json_report = output_file_path.replace(FILE_ENDING_QUALITY_FILTERED_FASTQ_GZ, FILE_ENDING_FASTP_JSON_REPORT)
    filepath_html_report = output_file_path.replace(FILE_ENDING_QUALITY_FILTERED_FASTQ_GZ, FILE_ENDING_FASTP_HTML_REPORT)

    #https://github.com/OpenGene/fastp/blob/59cc2f67414e74e99d42774e227b192a3d9bb63a/README.md#all-options
    command_fastp = [
        PROGRAM_PATH_FASTP, 
        "--thread", str(THREADS_DEFAULT),                    # Number of threads
        "--disable_adapter_trimming",
        "--qualified_quality_phred", "15",      #the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified.
        "--length_required" , "15",             #reads shorter than length_required will be discarded, default is 15. (int [=15])
        "--unqualified_percent_limit","40",     #how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% 
        "--n_base_limit", "5",                  #if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
        "--in1", input_file_path,               # Input R1 file
        "--out1", output_file_path,
        "--failed_out", filepath_failed_reads,
        "--json", filepath_json_report,
        "--html", filepath_html_report
    ]

    print_debug(f"Executing command: {' '.join(command_fastp)}")
    
    try:
        run_command(
            command_fastp, 
            description=f"fastp quality filter for {input_file_path}"
        )
        print_success(f"fastp quality filter for {input_file_path} complete")
    except subprocess.CalledProcessError as e:
        print_error(f"Failed to run fastp_quality_filter for {input_file_path}: {e}")


def fastp_quality_filter_for_species(species: str):

    print_species_execution(f"Running fastp quality filter for {species}")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_FILTERING,
        species
    ):
        print_skipping(f"Quality filtering is not enabled for {species} in the config. Skipping this step.")
        return

    reads_folder = get_folder_path_species_processed_adapter_removed(species)

    reads_files_list = get_files_in_folder_matching_pattern(reads_folder, f"*{FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ}")

    number_of_reads_files = len(reads_files_list)

    if number_of_reads_files == 0:
        print_warning(f"No adapter removed reads found for species {species}.")
        return
    
    print_debug(f"Found {len(reads_files_list)} adapter removed reads for species {species}.")
    print_debug(f"Adapter removed reads: {reads_files_list}")

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

        output_file_path = common_rrp.get_quality_filtered_path_for_adapter_removed_reads(species, read_file_path)
        execute_fastp_quality_filter(read_file_path, output_file_path)

    print_info(f"fastp quality filter for {species} complete")

def all_species_fastp_quality_filter():

    print_step_execution("Running fastp quality filter for all species")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.QUALITY_FILTERING
    ):
        print_skipping("Quality filtering is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        fastp_quality_filter_for_species(species)

    print_info("fastp quality filter for all species complete")