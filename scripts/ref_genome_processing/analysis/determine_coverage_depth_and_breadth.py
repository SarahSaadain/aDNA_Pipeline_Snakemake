import os
import subprocess
import collections
import pandas as pd

from multiprocessing import Pool, cpu_count

from common_aDNA_scripts import *

import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

def execute_samtools_detpth(input_file: str, coverage_output_file: str):

    print_info(f"Executing samtools depth for {input_file}")

    if not os.path.exists(input_file):
        raise Exception(f"Input file {input_file} does not exist!")

    if os.path.exists(coverage_output_file):
        print_skipping(f"Output file {coverage_output_file} already exists!")
        return

    command = f"{PROGRAM_PATH_SAMTOOLS} {PROGRAM_PATH_SAMTOOLS_DEPTH} -a {input_file} > {coverage_output_file}"
    print_debug(f"Samtools depth command: {command}")

    try:
        result = subprocess.run(command, shell=True,  check=True  )
        print_success(f"Samtools depth complete for {input_file}")
    except Exception as e:
        print_error(f"Failed to execute samtools depth: {e}")

def analyze_coverage_file(coverage_file, depth_breath_output_folder):

    # The `analyze_coverage_file()` function was originally processing the entire coverage file 
    # in-memory at once. However, this approach posed significant memory limitations when working 
    # with large coverage files (common in bioinformatics tasks). Processing large files as a whole 
    # would result in high memory usage and may cause the program to crash, particularly when 
    # dealing with files that contain millions of lines (such as BAM depth files).

    # ### Key Reasons for Switching to Chunked Processing:
    # 1. **Memory Efficiency:**
    #    - By processing the file in **chunks** rather than reading it all at once, we keep memory usage low. 
    #      This allows the program to handle much larger files, even if the system’s memory is limited.
    #    - The chunk size (set as `10^6` lines) can be adjusted based on system memory availability, offering flexibility.

    # 2. **Avoiding Memory Overload:**
    #    - If the entire file were read into memory at once (especially for large files), it could exhaust 
    #      available system RAM, leading to slowdowns or crashes. By processing chunks sequentially, only a 
    #      small portion of the file is in memory at any given time, reducing memory load.

    # 3. **Improved Performance on Large Datasets:**
    #    - Chunked processing allows us to process large datasets in manageable pieces. It helps avoid issues 
    #      with memory leaks or excessive garbage collection time (as only small portions of the file are processed at once).
    #    - It ensures that the program continues to run even for files that might otherwise exceed available memory.

    # 4. **Progress Tracking and Monitoring:**
    #    - When processing in chunks, progress can be tracked at meaningful intervals. 
    #      In this case, the program reports progress every time a chunk is processed, providing 
    #      real-time feedback on how much of the file has been completed.
    #    - This is particularly helpful in bioinformatics tasks, where large files may take a significant 
    #      amount of time to process.

    pid = os.getpid()  # Get current process ID for logging

    # Generate output file path based on input file name
    coverage_file_base_name = get_filename_from_path(coverage_file)
    analysis_file = coverage_file_base_name.replace(
        FILE_ENDING_SAMTOOLS_DEPTH_TSV,
        FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV
    )
    analysis_file_path = os.path.join(depth_breath_output_folder, analysis_file)

    print_info(f"Performing extended analysis for {coverage_file_base_name} with PID {pid}")

    # Skip processing if output already exists
    if os.path.exists(analysis_file_path):
        print_skipping(f"[PID {pid}] Analysis file {analysis_file_path} already exists.")
        return

    # Check that the input file exists
    if not os.path.exists(coverage_file):
        print_error(f"[PID {pid}] Coverage file {coverage_file} does not exist!")
        return

    # Count total lines in the file to estimate processing progress
    try:
        with open(coverage_file, 'r') as f:
            total_lines = sum(1 for _ in f)
    except Exception as e:
        print_error(f"[PID {pid}] Error counting lines in {coverage_file_base_name}: {e}")
        return
    
    print_info(f"[PID {pid}] Total lines in {coverage_file_base_name} to analyze: {total_lines:,}")

    # Dictionary to hold aggregated stats per scaffold
    summary_data = collections.defaultdict(lambda: {
        "depth_sum": 0,         # Sum of all depth values (used to calculate avg_depth)
        "max_depth": 0,         # Maximum depth seen
        "covered_bases": 0,     # Number of positions with depth > 0
        "total_bases": 0        # Total number of positions (regardless of depth)
    })

    lines_processed = 0                 # Running total of how many lines have been read
    chunk_size = 10**6                  # Number of rows per chunk (adjustable for tuning)

    # Dictionary to track progress milestones. Once a milestone is reached, it is marked as True
    # This prevents multiple printouts for the same milestone
    milestones = {25: False, 50: False, 75: False}

    try:
        # Process the file in chunks using pandas for memory efficiency
        for chunk in pd.read_csv(coverage_file, sep="\t", header=None,
                                 names=["scaffold", "position", "depth"],
                                 chunksize=chunk_size):

            # Group each chunk by scaffold and update cumulative statistics
            grouped = chunk.groupby("scaffold")
            for scaffold, group in grouped:
                summary_data[scaffold]["depth_sum"] += group["depth"].sum()
                summary_data[scaffold]["max_depth"] = max(
                    summary_data[scaffold]["max_depth"],
                    group["depth"].max()
                )
                summary_data[scaffold]["covered_bases"] += (group["depth"] > 0).sum()
                summary_data[scaffold]["total_bases"] += len(group)

            # Update line count and calculate percent progress
            lines_processed += len(chunk)
            percent_done = (lines_processed / total_lines) * 100

            # Log specific progress milestones (only once each)
            for milestone in milestones:
                if not milestones[milestone] and percent_done >= milestone:
                    print_info(f"[PID {pid}] Progress for {coverage_file_base_name}: {milestone}% completed")
                    milestones[milestone] = True

            # Regular progress printout (can be silenced if too verbose)
            print_debug(f"[PID {pid}] Progress for {coverage_file_base_name}: {lines_processed:,}/{total_lines:,} lines ({percent_done:.1f}%) processed")

    except Exception as e:
        print_error(f"[PID {pid}] Failed during processing of {coverage_file}: {e}")
        return

    # Convert cumulative summary data into a pandas DataFrame
    summary = pd.DataFrame.from_dict(summary_data, orient="index")
    summary.index.name = "scaffold"  # Label index for clarity

    # Calculate final metrics
    summary["avg_depth"] = summary["depth_sum"] / summary["total_bases"]
    summary["percent_covered"] = (summary["covered_bases"] / summary["total_bases"]) * 100

    # Reorder and filter columns to match expected output format
    summary = summary[["avg_depth", "max_depth", "covered_bases", "total_bases", "percent_covered"]]

    # Save result to CSV
    print_debug(f"[PID {pid}] Saving summary to {analysis_file_path} ...")
    summary.to_csv(analysis_file_path)

    print_info(f"[PID {pid}] Extended analysis complete for {coverage_file}")


def determine_coverage_depth_and_breath(species: str):
    """
    Orchestrates the coverage depth and breadth analysis for a single species.
    Executes samtools depth, performs extended analysis, and combines results.
    """
    print_species_execution(f"Processing coverage depth and breadth for species: {species}")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.COVERAGE_ANALYSIS,
        species
    ):
        print_skipping(f"Creating consensus sequence for species {species} is disabled in the config. Skipping this species.")
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

        print_info(f"Processing coverage depth and breadth for reference genome: {ref_genome_id}")

        # Step 1: Execute samtools depth for each BAM file
        execute_samtools_depth_for_species(species, ref_genome_id)

        # Step 2: Perform extended analysis on each coverage file
        perform_extended_analysis_for_species(species, ref_genome_id)

        # Step 3: Combine the individual analysis files
        combine_analysis_files(species, ref_genome_id)

    print_info(f"Coverage depth and breadth processing complete for species: {species}")

def execute_samtools_depth_for_species(species: str, reference_genome_id: str):
    """
    Finds all sorted BAM files for a species and executes samtools depth on each.
    """
    print_info(f"Executing samtools depth for all BAM files for species: {species}")

    mapped_folder = get_folder_path_species_processed_refgenome_mapped(species, reference_genome_id)
    list_of_bam_files = get_files_in_folder_matching_pattern(mapped_folder, f"*{FILE_ENDING_SORTED_BAM}")

    if len(list_of_bam_files) == 0:
        print_warning(f"No mapped BAM files found in {mapped_folder} for species {species}.")
        return

    print_debug(f"Found {len(list_of_bam_files)} BAM files for species {species}")
    print_debug(f"BAM files: {list_of_bam_files}")

    depth_breath_output_folder = get_folder_path_species_processed_refgenome_coverage(species, reference_genome_id)

    # Create a pool of worker processes to parallelize the execution.
    # Using cpu_count() to get the number of available CPU cores
    # use THREADS_DEFAULT to limit the number of processes as long as it is not higher than the number of cores
    # If THREADS_DEFAULT is higher than the number of cores, use the number of cores instead
    num_processes = min(THREADS_DEFAULT, cpu_count())
    
    with Pool(processes=num_processes) as pool:
        # Create a list of tuples, where each tuple contains the arguments
        # for the process_bam_file function.
        tasks = [(bam_file, depth_breath_output_folder) for bam_file in list_of_bam_files]
        # Use pool.starmap to apply the function to each set of arguments.
        # starmap unpacks the tuples and passes the elements as separate arguments.
        # process_bam_file is the function to be executed in parallel.
        # Each task is a tuple of (mapped_bam_file, depth_breath_output_folder)
        # The function will be called with the arguments from each tuple.
        # This allows for parallel processing of the BAM files.
        # starmap will block until all tasks are completed.
        pool.starmap(process_bam_file, tasks)

    print_info(f"Finished executing samtools depth for species {species}")

def process_bam_file(mapped_bam_file, depth_breath_output_folder):
    """
    Executes samtools depth on a single BAM file.  This helper function is for use with multiprocessing.
    """

    pid = os.getpid()

    mapped_bam_file_base_name = get_filename_from_path(mapped_bam_file)

    print_info(f"[PID {pid}] Processing BAM file {mapped_bam_file} with PID {pid}")

    coverage_file_name = mapped_bam_file_base_name.replace(FILE_ENDING_SORTED_BAM, FILE_ENDING_SAMTOOLS_DEPTH_TSV)
    coverage_output_file = os.path.join(depth_breath_output_folder, coverage_file_name)
    execute_samtools_detpth(mapped_bam_file, coverage_output_file)
    
    print_info(f"[PID {pid}] Finished samtools depth for {mapped_bam_file}") 

def perform_extended_analysis_for_species(species: str, reference_genome_id: str):
    """
    Finds all samtools depth output files for a species and performs extended analysis on each.
    """
    print_info(f"Performing extended analysis on depth files for species: {species}")

    samtools_depth_folder = get_folder_path_species_processed_refgenome_coverage(species, reference_genome_id)
    list_of_coverage_files = get_files_in_folder_matching_pattern(samtools_depth_folder, f"*{FILE_ENDING_SAMTOOLS_DEPTH_TSV}")

    if len(list_of_coverage_files) == 0:
        print_warning(f"No samtools depth files found in {samtools_depth_folder} for species {species}.")
        return

    print_debug(f"Found {len(list_of_coverage_files)} coverage files for species {species}")
    print_debug(f"Coverage files: {list_of_coverage_files}")

    target_folder = get_folder_path_species_results_refgenome_coverage(species, reference_genome_id)

    # Create a pool of worker processes to parallelize the execution.
    # Using cpu_count() to get the number of available CPU cores
    # use THREADS_DEFAULT to limit the number of processes as long as it is not higher than the number of cores
    # If THREADS_DEFAULT is higher than the number of cores, use the number of cores instead
    num_processes = min(THREADS_DEFAULT, cpu_count())

    with Pool(processes=num_processes) as pool:
        # Create a list of tuples, where each tuple contains the arguments
        # for the analyze_coverage_file function.
        tasks = [(file, target_folder) for file in list_of_coverage_files]
        # Use pool.starmap to apply the function to each set of arguments.
        # starmap unpacks the tuples and passes the elements as separate arguments.
        # analyze_coverage_file is the function to be executed in parallel.
        # Each task is a tuple of (coverage_file, depth_breath_output_folder)
        # The function will be called with the arguments from each tuple.
        # This allows for parallel processing of the coverage files.
        # starmap will block until all tasks are completed.
        pool.starmap(analyze_coverage_file, tasks)

    print_info(f"Finished performing extended analysis for species {species}")

def combine_analysis_files(species: str, reference_genome_id: str):
    
    print_info(f"Combining extended analysis files for species: {species}")

    # Folder paths
    individual_files_folder = get_folder_path_species_results_refgenome_coverage(species, reference_genome_id)
    analysis_folder = individual_files_folder  # same as above

    # Output file paths
    combined_file_path = os.path.join(analysis_folder, f"{species}{FILE_ENDING_COMBINED_COVERAGE_ANALYSIS_CSV}")
    combined_detailed_file_path = os.path.join(analysis_folder, f"{species}{FILE_ENDING_COMBINED_COVERAGE_ANALYSIS_DETAILED_CSV}")

    # Locate files
    individual_analysis_files = get_files_in_folder_matching_pattern(
        individual_files_folder,
        f"*{FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV}"
    )

    if not individual_analysis_files:
        print_warning(f"No individual analysis files found to combine for species {species} in {analysis_folder}.")
        return

    print_debug(f"Found {len(individual_analysis_files)} analysis files to combine for species {species}")
    print_debug(f"Analysis files: {individual_analysis_files}")

    combined_data = []       # For aggregated summary
    detailed_rows = []       # For raw, detailed per-scaffold data

    for analysis_file in individual_analysis_files:
        try:
            df_analysis = pd.read_csv(analysis_file)

            if df_analysis.empty:
                print_warning(f"Analysis file {analysis_file} is empty.")
                continue

            # Parse BAM filename from analysis filename
            original_bam_base = get_filename_from_path(analysis_file).replace(FILE_ENDING_EXTENDED_COVERAGE_ANALYSIS_CSV, "")
            original_bam_filename = original_bam_base + FILE_ENDING_SORTED_BAM

            # --- Aggregated stats ---
            total_bases_sum = df_analysis['total_bases'].sum()
            overall_avg_depth = (df_analysis['avg_depth'] * df_analysis['total_bases']).sum() / total_bases_sum if total_bases_sum > 0 else 0
            overall_max_depth = df_analysis['max_depth'].max()
            overall_covered_bases = df_analysis['covered_bases'].sum()
            overall_total_bases = total_bases_sum
            overall_percent_covered = (overall_covered_bases / overall_total_bases) * 100 if overall_total_bases > 0 else 0

            combined_data.append({
                "Filename": original_bam_filename,
                "OverallAvgDepth": overall_avg_depth,
                "OverallMaxDepth": overall_max_depth,
                "OverallCoveredBases": overall_covered_bases,
                "OverallTotalBases": overall_total_bases,
                "OverallPercentCovered": overall_percent_covered
            })

            # --- Append raw rows, tagged with filename ---
            df_analysis['Filename'] = original_bam_filename
            detailed_rows.append(df_analysis)

        except FileNotFoundError:
            print_error(f"Analysis file not found: {analysis_file}.")
        except pd.errors.EmptyDataError:
            print_warning(f"Analysis file is empty or malformed: {analysis_file}.")
        except Exception as e:
            print_error(f"An error occurred processing analysis file {analysis_file}: {e}")

    # --- Save aggregated summary ---
    if combined_data:
        df_combined = pd.DataFrame(combined_data)
        try:
            df_combined.to_csv(combined_file_path, index=False)
            print_success(f"Successfully created combined coverage analysis file: {combined_file_path}")
        except Exception as e:
            print_error(f"Error writing combined file: {e}")
    else:
        print_warning(f"No valid aggregated data to write for species {species}.")

    # --- Save detailed per-scaffold data ---
    if detailed_rows:
        try:
            df_detailed_combined = pd.concat(detailed_rows, ignore_index=True)
            df_detailed_combined.to_csv(combined_detailed_file_path, index=False)
            print_success(f"Successfully created detailed coverage file: {combined_detailed_file_path}")
        except Exception as e:
            print_error(f"Error writing detailed combined file: {e}")
    else:
        print_warning(f"No detailed data available to write for species {species}.")


def all_species_determine_coverage_depth_and_breath():

    print_step_execution("Determine coverage depth and breadth for all species")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.COVERAGE_ANALYSIS
    ):
        print_skipping("Determining coverage depth and breadth is not enabled in the config. Skipping this step.")
        return  

    for species in FOLDER_SPECIES: 
        determine_coverage_depth_and_breath(species)

    print_success("Finished determining coverage depth and breadth for all species")