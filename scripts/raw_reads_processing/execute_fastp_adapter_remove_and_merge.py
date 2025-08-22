import os
import subprocess
import raw_reads_processing.common_raw_reads_processing_helpers as common_rrp

from common_aDNA_scripts import *

def execute_fastp_paired_reads_remove_adapters_and_merge(input_file_path_r1: str, input_file_path_r2: str, output_file_path: str, adapter_sequence_r1:str, adapter_sequence_r2:str, threads:int = THREADS_DEFAULT):

    print_info(f"Removing adapters from {input_file_path_r1} and {input_file_path_r2} ...")

    if not os.path.exists(input_file_path_r1):
        raise Exception(f"Read file {input_file_path_r1[0]} does not exist!")

    if not os.path.exists(input_file_path_r2):
        raise Exception(f"Read file {input_file_path_r2} does not exist!")
    
    if os.path.exists(output_file_path):
        print_info(f"Output file {output_file_path} already exists! Skipping!")
        return
    
    filepath_merge_failed_passed_r1 = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, FILE_ENDING_MERGE_FAILED_PASSED_R1_FASTQ_GZ)
    filepath_merge_failed_passed_r2 = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, FILE_ENDING_MERGE_FAILED_PASSED_R2_FASTQ_GZ)
    filepath_merge_failed_not_passed_r1 = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, FILE_ENDING_MERGE_FAILED_NOT_PASSED_R1_FASTQ_GZ)
    filepath_merge_failed_not_passed_r2 = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, FILE_ENDING_MERGE_FAILED_NOT_PASSED_R2_FASTQ_GZ)
    filepath_merge_json_report = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, FILE_ENDING_FASTP_JSON_REPORT)
    filepath_merge_html_report = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, FILE_ENDING_FASTP_HTML_REPORT)

    #https://github.com/OpenGene/fastp/blob/59cc2f67414e74e99d42774e227b192a3d9bb63a/README.md#all-options
    command_fastp = [
        PROGRAM_PATH_FASTP,
        "--adapter_sequence", adapter_sequence_r1,  # Adapter for R1
        "--adapter_sequence_r2", adapter_sequence_r2,  # Adapter for R2
        "--out1", filepath_merge_failed_passed_r1,
        "--out2", filepath_merge_failed_passed_r2,
        "--unpaired1", filepath_merge_failed_not_passed_r1,
        "--unpaired2", filepath_merge_failed_not_passed_r2,
        "--merged_out", output_file_path,  # Output file for R1
        "--in1", input_file_path_r1,        # Input R1 file
        "--in2", input_file_path_r2,         # Input R2 file
        "--json", filepath_merge_json_report,
        "--html", filepath_merge_html_report,
        "--merge",
        "--thread", str(threads),               # Number of threads

        # length filtering options
        "--length_required" , "15",             #reads shorter than length_required will be discarded, default is 15. (int [=15])
        
        # poly
        "--trim_poly_x", "5",

        # quality filtering options
        "--qualified_quality_phred", "5",       #the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
        "--unqualified_percent_limit", "40",    #how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])"
        "--n_base_limit", "5"                  #if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
       
    ]

    print_debug(f"Executing command: {' '.join(command_fastp)}")
    
    try:
        subprocess.run(command_fastp, check=True)
        print_success(f"Adapters removed from {input_file_path_r1} and {input_file_path_r2}.")
    except subprocess.CalledProcessError as e:
        raise Exception(f"Removed adapters error for {input_file_path_r1} and {input_file_path_r2} : {e}")


def execute_fastp_single_reads_remove_adapters(input_file_path: str, output_file_path: str, adapter_sequence: str, threads: int = THREADS_DEFAULT):
    
    print_info(f"Removing adapters from {input_file_path} ...")
    
    if not os.path.exists(input_file_path):
        raise Exception(f"Read file {input_file_path} does not exist!")
    
    if os.path.exists(output_file_path):
        print_skipping(f"Output file {output_file_path} already exists!")
        return
    
    filepath_json_report = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, "_report.json")
    filepath_html_report = output_file_path.replace(FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, "_report.html")
    
    command_fastp = [
        PROGRAM_PATH_FASTP,
        "--adapter_sequence", adapter_sequence,  # Adapter sequence for single-end reads
        "-i", input_file_path,  # Input file
        "-o", output_file_path,  # Output file
        "--json", filepath_json_report,
        "--html", filepath_html_report,
        "--thread", str(threads),  # Number of threads

        # length filtering options
        "--length_required", "15",  # Reads shorter than this will be discarded

        # poly
        "--trim_poly_x", "5",

        # quality filtering options
        "--qualified_quality_phred", "5",  # Minimum quality threshold
        "--unqualified_percent_limit", "40",  # Max percentage of unqualified bases
        "--n_base_limit", "5"  # Max number of N bases allowed
    ]

    print_debug(f"Executing command: {' '.join(command_fastp)}")
    
    try:
        run_command(
            command_fastp, 
            description=f"fastp dapters removal for {input_file_path}"
        )

        print_success(f"Adapters removed from {input_file_path}.")
    except subprocess.CalledProcessError as e:
        raise Exception(f"Adapter removal error for {input_file_path}: {e}")


def all_species_fastp_adapter_remove_and_merge():

    print_step_execution("Running adaper removal and merge for all species")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.ADAPTER_REMOVAL
    ):
        print_skipping("Adapter removal and merge is not enabled in the config. Skipping this step.")
        return
    
    for species in FOLDER_SPECIES:
        adapter_remove_for_species(species)

    print_info("Adapter removal and merge for all species completed successfully.")

def adapter_remove_for_species(species: str):
    print_species_execution(f"Running adapter removal for species {species}")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.ADAPTER_REMOVAL,
        species
    ):
        print_skipping(f"Adapter removal is not enabled for {species} in the config. Skipping this step.")
        return

    # Get lists of R1 and R2 read files
    raw_reads_folder = get_folder_path_species_raw_reads(species)
    list_of_r1_read_files = get_files_in_folder_matching_pattern(raw_reads_folder, FILE_PATTERN_R1_FASTQ_GZ)
    list_of_r2_read_files = get_files_in_folder_matching_pattern(raw_reads_folder, FILE_PATTERN_R2_FASTQ_GZ)

    if not list_of_r1_read_files and not list_of_r2_read_files:
        print_warning(f"No raw reads found for species {species}.")
        return

    # Create a set of base filenames to identify paired reads
    paired_reads = set()
    single_reads = set(list_of_r1_read_files)  # Assume all R1 files are single initially

    for r1 in list_of_r1_read_files:
        r2 = r1.replace("_R1_", "_R2_")  # Generate expected R2 filename
        if r2 in list_of_r2_read_files:
            paired_reads.add((r1, r2))
            single_reads.discard(r1)  # Remove from single-end list if paired

    #get adapter sequence
    adapter_sequence_r1, adapter_sequence_r2 = get_adapter_sequence(species)

    # Process paired reads
    if paired_reads:

        number_of_paired_reads = len(paired_reads)
        count_current = 0

        print_debug(f"Found {number_of_paired_reads} paired-end reads for species {species}.")
        print_debug(f"Paired-end reads: {paired_reads}")

        print_info("Processing paired-end reads.")
        try:
            for r1, r2 in paired_reads:

                count_current += 1

                print_info(f"[{count_current}/{number_of_paired_reads}] Processing paired reads: {get_filename_from_path(r1)}, {get_filename_from_path(r2)}")  

                adapter_removed_read_file = common_rrp.get_adapter_removed_path_for_paired_raw_reads(species, [r1, r2])
                
                individual = common_rrp.get_individual_from_file(r1)

                is_ref_genome_read_file_exists = common_rrp.is_species_individual_reads_file_exists(
                    species, 
                    individual
                )

                if is_ref_genome_read_file_exists:
                    print_skipping(f"Individual {individual} already prepared for reference genome processing!")
                    continue
                
                execute_fastp_paired_reads_remove_adapters_and_merge(r1, r2, adapter_removed_read_file, adapter_sequence_r1, adapter_sequence_r2)
       
        except Exception as e:
            print_error(f"Error processing paired-end reads for species {species}: {e}")
            return
    
    # Process single-end reads
    if single_reads:

        number_of_single_reads = len(single_reads)
        count_current = 0

        print_debug(f"Found {number_of_single_reads} single-end reads for species {species}.")
        print_debug(f"Single-end reads: {single_reads}")
        
        print_info("Processing single-end reads.")
        try:
            for read_file_path in single_reads:

                count_current += 1

                print_info(f"[{count_current}/{number_of_single_reads}] Processing single read: {get_filename_from_path(read_file_path)}")

                adapter_removed_read_file = common_rrp.get_adapter_removed_path_for_paired_raw_reads(species, [read_file_path])
                
                individual = common_rrp.get_individual_from_file(read_file_path)

                is_ref_genome_read_file_exists = common_rrp.is_species_individual_reads_file_exists(
                    species, 
                    individual
                )

                if is_ref_genome_read_file_exists:
                    print_skipping(f"Individual {individual} already prepared for reference genome processing!")
                    continue
                
                execute_fastp_single_reads_remove_adapters(read_file_path, adapter_removed_read_file, adapter_sequence_r1)
        except Exception as e:
            print_error(f"Error processing single-end reads for species {species}: {e}")
            return
    
    print_info(f"Adapter removal for species {species} completed successfully.")
