import os
import pandas as pd
import gzip
import raw_reads_processing.common_raw_reads_processing_helpers as common_rrp

from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio import SeqIO
from common_aDNA_scripts import *

def get_read_length_distribution(fastq_file: str) -> Counter:
    read_lengths = Counter()

    with gzip.open(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            read_lengths[len(record.seq)] += 1  # Count occurrences of each read length

    return read_lengths
    
def get_file_name_read_length_distribution(species: str) -> str:
    return f"{species}{FILE_ENDING_READ_LENGTH_DISTRIBUTION_TSV}"

def _process_single_read_length_file(raw_read_path: str, species: str) -> str | None:

    try:
        pid = os.getpid() # Get current process ID for logging
        print_info(f"[PID {pid}] Processing read file {raw_read_path} for read length distribution")

        # Define the directory where this individual file's results will be saved.
        processed_dir = get_folder_path_species_processed_qc_read_length_distribution(species)

        # Extract metadata from the filename to create a unique temporary file name.
        reads_file_id = get_filename_from_path_without_extension(raw_read_path)
        
        temp_file_name = f"{reads_file_id}{FILE_ENDING_READ_LENGTH_DISTRIBUTION_TSV}"
        temp_file_path = os.path.join(processed_dir, temp_file_name)

        # Check if the individual result file already exists to avoid redundant processing.
        if os.path.exists(temp_file_path):
            print_skipping(f"[PID {pid}] Individual result file already exists: {temp_file_path}.")
            return temp_file_path # Return existing path if already processed

        # Get read length distribution for raw reads
        raw_distribution = get_read_length_distribution(raw_read_path)
        
        # Get read length distribution for adapter-removed reads
        # get_adapter_removed_path_for_paired_raw_reads expects a list, so wrap raw_read_path
        adapter_removed_file = common_rrp.get_adapter_removed_path_for_paired_raw_reads(species, [raw_read_path])
        print_info(f"[PID {pid}] Processing adapter removed file {adapter_removed_file}")
        adapter_removed_distribution = get_read_length_distribution(adapter_removed_file)

        # Get read length distribution for quality-filtered reads
        quality_filtered_file = common_rrp.get_quality_filtered_path_for_adapter_removed_reads(species, adapter_removed_file)
        print_info(f"[PID {pid}] Processing quality filtered file {quality_filtered_file}")
        quality_filtered_distribution = get_read_length_distribution(quality_filtered_file)

        # Get read length distribution for deduplicated reads
        duplicates_removed_file = common_rrp.get_deduplication_path_for_quality_filtered_reads(species, quality_filtered_file)
        print_info(f"[PID {pid}] Processing duplicates removed file {duplicates_removed_file}")
        duplicates_removed_distribution = get_read_length_distribution(duplicates_removed_file)

        # Extract metadata (individual, protocol) from the reads_file_id.
        # Assumes a filename format like "individual_protocol_R1.fastq.gz"
        parts = reads_file_id.split("_")
        individual = parts[0] if len(parts) > 0 else "N/A"
        protocol = parts[1] if len(parts) > 1 else "N/A"

        # Convert Counter objects to pandas DataFrames
        df_raw = pd.DataFrame(raw_distribution.items(), columns=["read_length", "read_count_raw"])
        df_adapter = pd.DataFrame(adapter_removed_distribution.items(), columns=["read_length", "read_count_adapter_removed"])
        df_quality = pd.DataFrame(quality_filtered_distribution.items(), columns=["read_length", "read_count_quality_filtered"])
        df_dedup = pd.DataFrame(duplicates_removed_distribution.items(), columns=["read_length", "read_count_duplicates_removed"])

        # Merge all DataFrames on 'read_length' using an outer join.
        # This ensures all read lengths present in any stage are included, filling missing counts with 0.
        df = df_raw.merge(df_adapter, on="read_length", how="outer") \
                .merge(df_quality, on="read_length", how="outer") \
                .merge(df_dedup, on="read_length", how="outer")

        # Fill any NaN values (from outer merge where a length was not present in all stages) with 0
        # and convert all count columns to integer type.
        df = df.fillna(0).astype(int)

        # Add metadata columns to the DataFrame, inserting them at the beginning.
        df.insert(0, "protocol", protocol)
        df.insert(0, "individual", individual)
        df.insert(0, "reads_file", reads_file_id)
        
        # Save this single-row DataFrame to its dedicated temporary TSV file.
        df.to_csv(temp_file_path, sep="\t", index=False)
        
        print_info(f"[PID {pid}] Saved individual result file: {temp_file_path}")
    except Exception as e:
        print_error(f"[PID {pid}] Error processing {raw_read_path} in parallel: {e}")


def determine_read_length_distribution(species: str):

    print_species_execution(f"Determine read length distribution for {species}")

    # Get a list of all raw read files (R1) for the specified species.
    raw_reads = get_files_in_folder_matching_pattern(get_folder_path_species_raw_reads(species), FILE_PATTERN_R1_FASTQ_GZ)

    # If no raw reads are found, log a warning and exit the function.
    if len(raw_reads) == 0:
        print_warning(f"No raw reads found for species {species}.")
        return
    
    print_debug(f"Found {len(raw_reads)} raw reads for species {species}.")
    print_debug(f"Raw reads: {raw_reads}")

    # we get the processed directory for the species in order to create the folder before
    # starting the parallel processing
    # This is important to ensure that the directory exists before any worker tries to write to it.
    processed_dir = get_folder_path_species_processed_qc_read_length_distribution(species)
    
    # Determine the number of processes to use for parallel execution.
    # It takes the minimum of a predefined default (THREADS_DEFAULT) and the
    # actual number of CPU cores available on the system, ensuring efficient
    # resource utilization without overloading the system.
    num_processes = min(THREADS_DEFAULT, cpu_count())
    print_debug(f"Using {num_processes} processes for parallel execution.")

    # Initialize a multiprocessing Pool. This creates a pool of worker processes
    # that can execute tasks concurrently. The 'processes' argument specifies
    # the maximum number of worker processes to use.
    with Pool(processes=num_processes) as pool:
        # Prepare the arguments for each call to the _process_single_read_length_file helper function.
        # Each tuple (raw_read, species) represents one set of arguments for a single task.
        args_for_pool = [(raw_read, species) for raw_read in raw_reads]

        # Use pool.starmap to apply the _process_single_read_length_file function to each
        # set of arguments in args_for_pool. starmap is suitable when the target
        # function expects multiple arguments, which are provided as a tuple.
        # This call blocks until all tasks in the pool have completed, and returns
        # a list of results (paths to saved files or None if an error occurred).
        pool.starmap(_process_single_read_length_file, args_for_pool)


def combine_read_length_distributions(species: str):
        
        print_info(f"Combining read length distribution for {species}")

        # Get a list of all temporary files created during the parallel processing.
        all_temp_files = get_files_in_folder_matching_pattern(
            get_folder_path_species_processed_qc_read_length_distribution(species), 
            f"*{FILE_ENDING_READ_LENGTH_DISTRIBUTION_TSV}"
        )
    
        # If no individual result files were successfully created, log a warning and exit.
        if not all_temp_files:
            print_warning(f"No read length distribution data collected for species {species}. This might indicate errors during parallel processing or no raw reads were found.")
            return
        
        print_debug(f"Found {len(all_temp_files)} individual result files for species {species}.")
        print_debug(f"Individual result files: {all_temp_files}")

        print_info(f"Combining {len(all_temp_files)} individual result files into a single output.")
        # Initialize an empty DataFrame to collect all results.
        # Columns are explicitly defined to ensure order and presence even if no data is collected.
        combined_df = pd.DataFrame(columns=[
            "reads_file", "individual", "protocol", "read_length",
            "read_count_raw", "read_count_adapter_removed",
            "read_count_quality_filtered", "read_count_duplicates_removed"
        ])

        # Define the final output file path for the combined results.
        output_file_path = os.path.join(get_folder_path_species_results_qc_read_length_distribution(species),  get_file_name_read_length_distribution(species))

        # Iterate through each successfully created temporary file.
        for temp_file in all_temp_files:
            try:
                # Read the content of each temporary file into a small DataFrame.
                # These files are expected to contain the read length distribution for one sample.
                df_part = pd.read_csv(temp_file, sep="\t")
                # Concatenate the individual DataFrame to the main combined DataFrame.
                # ignore_index=True ensures that the new rows are appended with a continuous index.
                combined_df = pd.concat([combined_df, df_part], ignore_index=True)
            except Exception as e:
                # Log a warning if there's an error reading any individual file,
                # but continue processing other files.
                print_warning(f"Error reading individual result file {temp_file} during combination: {e}")

        # Save the final combined DataFrame to the specified output file in TSV format.
        # index=False prevents pandas from writing the DataFrame index as a column.
        combined_df.to_csv(output_file_path, sep="\t", index=False)
        print_info(f"Successfully combined results and saved to: {output_file_path}") 


def all_species_determine_read_length_distribution():
    print_step_execution("Determine reads processing result for all species")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.READ_LENGTH_DISTRIBUTION_ANALYSIS
    ):
        print_skipping(f"Reads length distribution analysis is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 

        if not config.is_process_step_enabled(
            PipelineStages.RAW_READS_PROCESSING, 
            RawReadsProcessingSteps.READ_LENGTH_DISTRIBUTION_ANALYSIS,
            species
        ):
            print_skipping(f"Reads length distribution analysis for {species} is not enabled in the config. Skipping this step.")
            continue

        determine_read_length_distribution(species)
        combine_read_length_distributions(species)
        print_success(f"Finished determining reads processing result for {species}")

    print_success("Finished determining reads processing result for all species")