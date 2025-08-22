import os
import re
import subprocess
import pandas as pd
import raw_reads_processing.common_raw_reads_processing_helpers as common_rrp

from multiprocessing import Pool, cpu_count
from common_aDNA_scripts import *

def execute_seqkit_stats_count_reads(input_file: str, thread:int = THREADS_DEFAULT) -> int:

    print_info(f"Executing seqkit stats for {input_file}")

    if not os.path.exists(input_file):
        raise Exception(f"Input file {input_file} does not exist!")

    try:
        command = [PROGRAM_PATH_SEQKIT, PROGRAM_PATH_SEQKIT_STATS, "--threads", str(thread), input_file]

        print_debug(f"Executing command: {' '.join(command)}")

        output = run_command(
            command, 
            description=f"Executing seqkit stats for {input_file}",
            throw_error=True
        )
    
    except Exception as e:
        print_error(f"Failed to execute seqkit stats: {e}")
        return -1

    try:

        # the result looks like this:
        # seqkit stats -j 10 /mnt/data2/sarah/aDNA/Mmus/raw/reads/326862_S37_R1_001.fastq.gz
        # file                                                             format  type   num_seqs        sum_len  min_len  avg_len  max_len
        # /mnt/data2/sarah/aDNA/Mmus/raw/reads/326862_S37_R1_001.fastq.gz  FASTQ   DNA   7,729,611  1,159,441,650      150      150      150

        # Read output into pandas DataFrame so that we can get the data we need

        # Print raw output to debug the format

        # Define the regular expression pattern to capture columns
        pattern = r"([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([0-9,]+)\s+([0-9,]+)\s+([0-9,]+)\s+([0-9,\.]+)\s+([0-9,]+)"

        # Use regex to extract all lines matching the pattern
        matches = re.findall(pattern, output)

        # If no matches were found, raise an error
        if not matches:
            raise ValueError("No valid output found to parse.")

        # Convert matches to a DataFrame
        df = pd.DataFrame(matches, columns=["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"])

        # Convert numeric columns to correct types
        df[["num_seqs", "sum_len", "min_len", "max_len"]] = df[["num_seqs", "sum_len", "min_len", "max_len"]].apply(lambda x: x.str.replace(',', '').astype(int))

        # Extract number of sequences from the DataFrame
        sequences_count = df["num_seqs"].iloc[0]

        print_debug(f"Seqkit stats complete for {input_file}: {sequences_count} reads")

        return sequences_count
    except Exception as e:
        print_error(f"Failed to process seqkit stats output: {e}")
        return -1
    
def get_file_name_reads_processing(species: str) -> str:
    return f"{species}{FILE_ENDING_READS_PROCESSING_RESULT_TSV}"

# --- Helper function for parallel processing ---
def _process_single_read_file(raw_read_path: str, species: str):
    try:

        pid = os.getpid()  # Get current process ID for logging

        print_info(f"[PID {pid}] Processing single read file: {raw_read_path}")

         # Create a temporary directory to store individual results from parallel processes
        processed_dir = get_folder_path_species_processed_qc_reads_processing(species)

        # Extract reads_id, individual, and protocol from the filename
        reads_id = get_filename_from_path_without_extension(raw_read_path)
        
        temp_file_name = f"{reads_id}{FILE_ENDING_READS_PROCESSING_RESULT_TSV}" # Using .tsv for tab-separated
        temp_file_path = os.path.join(processed_dir, temp_file_name)

        if os.path.exists(temp_file_path):
            print_skipping(f"[PID {pid}] File already exists: {temp_file_path}.")
            return

        # Count raw reads
        raw_count = execute_seqkit_stats_count_reads(raw_read_path, thread=1)

        # Determine paths and count reads after adapter removal
        # Note: get_adapter_removed_path_for_paired_raw_reads expects a list,
        # so we wrap raw_read_path in a list for compatibility with its placeholder.
        print_info(f"[PID {pid}] Processing adapter removed file for {reads_id}")
        adapter_removed_file = common_rrp.get_adapter_removed_path_for_paired_raw_reads(species, [raw_read_path])
        adapter_removed_count = execute_seqkit_stats_count_reads(adapter_removed_file, thread=1)
        print_info(f"[PID {pid}] Adapter removed count: {adapter_removed_count}")

        # Determine paths and count reads after quality filtering
        print_info(f"[PID {pid}] Processing quality filtered file for {reads_id}")
        quality_filtered_file = common_rrp.get_quality_filtered_path_for_adapter_removed_reads(species, adapter_removed_file)
        quality_filtered_count = execute_seqkit_stats_count_reads(quality_filtered_file, thread=1)
        print_info(f"[PID {pid}] Quality filtered count: {quality_filtered_count}")

        # Determine paths and count reads after deduplication
        print_info(f"[PID {pid}] Processing duplicates removed file for {reads_id}")
        duplicates_removed_file = common_rrp.get_deduplication_path_for_quality_filtered_reads(species, quality_filtered_file)
        duplicates_removed_count = execute_seqkit_stats_count_reads(duplicates_removed_file, thread=1)
        print_info(f"[PID {pid}] Duplicates removed count: {duplicates_removed_count}")

        # Assuming filename format like "individual_protocol_R1.fastq.gz"
        parts = reads_id.split("_")
        individual = parts[0] if len(parts) > 0 else "N/A"
        protocol = parts[1] if len(parts) > 1 else "N/A"

        # Create a DataFrame for the current row
        new_row_df = pd.DataFrame({
            "reads_file": [reads_id],
            "individual": [individual],
            "protocol": [protocol],
            "raw_count": [raw_count],
            "adapter_removed_count": [adapter_removed_count],
            "quality_filtered_count": [quality_filtered_count],
            "duplicates_removed_count": [duplicates_removed_count]
        })

        # Save this single row to a temporary file
        
        new_row_df.to_csv(temp_file_path, sep="\t", index=False)

        print_info(f"[PID {pid}] Saved file: {temp_file_path}")
    except Exception as e:
        print_error(f"[PID {pid}] Error processing {raw_read_path}: {e}")

def determine_reads_processing_result(species: str):

    print_species_execution(f"Determine reads processing result for {species}")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.READS_PROCESSING_ANALYSIS,
        species
    ):
        print_skipping(f"Reads processing analysis is not enabled for {species} in the config. Skipping this step.")
        return

    # Get a list of raw read files for the species
    raw_reads = get_files_in_folder_matching_pattern(
        get_folder_path_species_raw_reads(species),
        FILE_PATTERN_R1_FASTQ_GZ
    )

    if len(raw_reads) == 0:
        print_warning(f"No raw reads found for species {species}.")
        return

    print_debug(f"Found {len(raw_reads)} raw reads for species {species}.")
    print_debug(f"Raw reads: {raw_reads}")

    # Determine the number of processes to use for parallel execution.
    # It takes the minimum of a predefined default (THREADS_DEFAULT) and the
    # actual number of CPU cores available on the system, ensuring efficient
    # resource utilization without overloading the system.
    num_processes = min(THREADS_DEFAULT, cpu_count())
    print_debug(f"Using {num_processes} processes for parallel execution.")

    # we get the processed directory for the species in order to create the folder before
    # starting the parallel processing
    # This is important to ensure that the directory exists before any worker tries to write to it.
    processed_dir = get_folder_path_species_processed_qc_reads_processing(species)

    # Initialize a multiprocessing Pool. This creates a pool of worker processes
    # that can execute tasks concurrently. The 'processes' argument specifies
    # the maximum number of worker processes to use.
    with Pool(processes=num_processes) as pool:
        # Prepare the arguments for each call to the _process_single_read_file helper function.
        # Each tuple (raw_read, species) represents one set of arguments for a single task.
        args_for_pool = [(raw_read, species) for raw_read in raw_reads]

        # Use pool.starmap to apply the _process_single_read_file function to each
        # set of arguments in args_for_pool. starmap is suitable when the target
        # function expects multiple arguments, which are provided as a tuple.
        # This call blocks until all tasks in the pool have completed.
        pool.starmap(_process_single_read_file, args_for_pool)

def combine_reads_processing_results(species: str):

    print_info(f"Combining reads processing results for {species}")

     # Get a list of raw read files for the species
    analysis_files = get_files_in_folder_matching_pattern(
        get_folder_path_species_processed_qc_reads_processing(species),
        f"*{FILE_ENDING_READS_PROCESSING_RESULT_TSV}"
    )

    if not analysis_files:
        print_warning(f"No analysis files found for species {species}.")
        return
    
    print_debug(f"Found {len(analysis_files)} analysis files for species {species}.")
    print_debug(f"Analysis files: {analysis_files}")

    # Initialize an empty DataFrame to collect all results
    # Define columns explicitly to ensure order and presence even if no data is collected
    combined_df = pd.DataFrame(columns=[
        "reads_file", "individual", "protocol", "raw_count",
        "adapter_removed_count", "quality_filtered_count", "duplicates_removed_count"
    ])

        # Define the final output file path
    output_file_path = os.path.join(
        get_folder_path_species_results_qc_reads_processing(species),
        get_file_name_reads_processing(species)
    )

    # Combine all individual temporary result files into a single DataFrame
    print_info(f"Combining {len(analysis_files)} individual result files...")
    for temp_file in analysis_files:
        try:
            # Read each temporary file (which contains a single row)
            df_part = pd.read_csv(temp_file, sep="\t")
            # Concatenate it to the main combined DataFrame
            combined_df = pd.concat([combined_df, df_part], ignore_index=True)
        except Exception as e:
            print_warning(f"Error reading temporary file {temp_file} during combination: {e}")

    # Save the final combined DataFrame to the specified output file
    combined_df.to_csv(output_file_path, sep="\t", index=False)
    print_info(f"Successfully combined results and saved to: {output_file_path}")

def all_species_determine_determine_reads_processing_result():
    print_step_execution("Determine reads processing result for all species")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.READS_PROCESSING_ANALYSIS
    ):
        print_skipping(f"Reads processing analysis is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 

        if not config.is_process_step_enabled(
            PipelineStages.RAW_READS_PROCESSING, 
            RawReadsProcessingSteps.READS_PROCESSING_ANALYSIS,
            species
        ):
            print_skipping(f"Reads processing analysis for {species} is not enabled in the config. Skipping this step.")
            continue

        determine_reads_processing_result(species)
        combine_reads_processing_results(species)
        print_info(f"Finished determining reads processing result for {species}")

    print_success("Finished determining reads processing result for all species")