import os
import subprocess
import pandas as pd
from collections import defaultdict

from common_aDNA_scripts import *


def run_kraken_on_file(species: str, fastq_file_path: str,  Kraken_report_tsv: str, threads: int = THREADS_DEFAULT):
   
    print_info(f"Running Kraken on file: {os.path.basename(fastq_file_path)}")

    # Check if output files already exist
    if os.path.exists(Kraken_report_tsv):
        print_skipping(f"Output files {Kraken_report_tsv} already exist.")
        return
    
    # get the Kraken database path from the config
    Kraken_db = config.get_pipeline_settings(RawReadsProcessingSteps.CONTAMINATION_CHECK).get(ContaminationCheckSettings.KRAKEN_DB.value)

    # Ensure Kraken2 database path is set and exists
    if not Kraken_db:
        print_error("Kraken2 database path is not set. Please check your pipeline configuration.")
        return
    
    if not os.path.exists(Kraken_db):
        print_error(f"Kraken2 database path does not exist: {Kraken_db}")
        return

    # Construct the kraken2 command
    # Using --gzip-compressed assumes the input is .gz
    kraken2_command = [
        PROGRAM_PATH_KRAKEN,
        "--db", Kraken_db,
        "--threads", str(threads),
        "--gzip-compressed",
        "--output", Kraken_report_tsv,
        fastq_file_path # Input file
    ]

    print_debug(f"Kraken2 command: {' '.join(kraken2_command)}")

    # Execute the command
    try:
        # Using capture_output=True and text=True to get stdout/stderr in case of errors
        result = subprocess.run(kraken2_command, check=True, capture_output=True, text=True)
        print_success(f"Kraken2 analysis complete for {get_filename_from_path(fastq_file_path)}")
        # Optionally print stdout/stderr for debugging
        print_debug("Kraken2 stdout:\n" + result.stdout)
        print_debug("Kraken2 stderr:\n" + result.stderr)
    except subprocess.CalledProcessError as e:
        print_error(f"Kraken2 failed for {get_filename_from_path(fastq_file_path)} with error: {e.returncode}")
        print_error("Kraken2 stdout:\n" + e.stdout)
        print_error("Kraken2 stderr:\n" + e.stderr)
    except FileNotFoundError:
         print_error("Kraken2 command not found. Make sure kraken2 is installed and in your PATH.")
    except Exception as e:
        print_error(f"An unexpected error occurred while running Kraken2 on {get_filename_from_path(fastq_file_path)}: {e}")

def create_kraken_top5_analysis(report_file_path: str, output_file_path: str):

    print_info(f"Analyzing Kraken2 report: {get_filename_from_path(report_file_path)}")

    # Check if output file already exists
    if os.path.exists(output_file_path):
        print_skipping(f"Analysis output file {get_filename_from_path(output_file_path)} already exists.")
        return

    # Construct the shell command pipeline
    # awk '$1 == "C" {print $3}': Filters lines starting with 'C' and prints the 3rd column (species name)
    # sort: Sorts the species names alphabetically
    # uniq -c: Counts occurrences of each unique species name
    # sort -nr: Sorts the counts in reverse numerical order
    # head -5: Takes the top 5 entries
    analysis_command = f"awk '$1 == \"C\" {{print $3}}' {report_file_path} | sort | uniq -c | sort -nr | head -5 > {output_file_path}"

    print_debug(f"Analysis command: {analysis_command}")

    # Execute the command
    try:
        # Use shell=True because we are using a pipeline with pipes (|) and redirection (>)
        subprocess.run(analysis_command, shell=True, check=True, capture_output=True, text=True)
        print_success(f"Analysis complete. Top 5 species written to {get_filename_from_path(output_file_path)}")
        # Optionally print stdout/stderr for debugging
        # print_debug("Analysis stdout:\n" + result.stdout)
        # print_debug("Analysis stderr:\n" + result.stderr)
    except subprocess.CalledProcessError as e:
        print_error(f"Analysis failed for {get_filename_from_path(report_file_path)} with error: {e.returncode}")
        print_error("Analysis stdout:\n" + e.stdout)
        print_error("Analysis stderr:\n" + e.stderr)
    except Exception as e:
        print_error(f"An unexpected error occurred during analysis of {get_filename_from_path(report_file_path)}: {e}")

def combine_kraken2_top5_analysis(species: str):
    # Print which species is being analyzed
    print(f"Combining Kraken2 top 5 analysis for species: {species}")
    
    # Get the folder where the Kraken2 result files are stored
    results_folder = get_folder_path_species_results_qc_kraken(species)

    # Get a list of all files in the folder that match the Kraken2 TSV pattern
    kraken_reports = get_files_in_folder_matching_pattern(
        results_folder,
        f"*{FILE_ENDING_KRAKEN_TOP5_ANALYSIS_TSV}"
    )

    # If no Kraken2 reports are found, print a warning and return
    if not kraken_reports:
        print_warning(f"No Kraken2 top 5 analysis found for species {species} in {results_folder}. Skipping analysis.")
        return
    
    print_debug(f"Found {len(kraken_reports)} Kraken2 top 5 analysis for species {species}")
    print_debug(f"Files found: {[get_filename_from_path(f) for f in kraken_reports]}")

    # Dictionary to store the count of each species per file
    all_counts = defaultdict(dict)  # Format: {filename: {species_id: count}}

    # Dictionary to store total count of each species across all files
    # This will be used to sort species by total count
    species_total_counts = defaultdict(int)  # Format: {species_id: total_count}

    # Loop through each Kraken2 file
    for filepath in kraken_reports:
        filename = get_filename_from_path(filepath)  # Get the filename from the full path

        # Open and read the file line by line
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) != 2:
                    continue  # Skip lines that don’t have exactly two values
                count, species_id = map(int, parts)  # Convert both values to integers

                # Store count for the species in the current file
                all_counts[filename][species_id] = count

                # Add count to total across all files
                species_total_counts[species_id] += count

    # If no species were found, print a warning and return
    if not species_total_counts:
        print_warning(f"No contaminants found in Kraken2 top 5 analysis for species {species}.")
        return
    
    print_debug(f"Total species found: {len(species_total_counts)}")
    print_debug(f"Species total counts: {species_total_counts}")

    # Sort species IDs based on their total count (descending)
    sorted_species_ids = sorted(species_total_counts, key=species_total_counts.get, reverse=True)

    # Create a list of rows for the final DataFrame
    rows = []
    for filename, species_counts in all_counts.items():
        
        parts = filename.split('_')
        individuum = parts[0]   # e.g., "Bger3"
        protocol = parts[1]     # e.g., "D"

        row = {
            'filename': filename,
            'individuum': individuum,
            'protocol': protocol
        }

        for species_id in sorted_species_ids:
            # Fill in count if it exists, otherwise put 0
            row[species_id] = species_counts.get(species_id, 0)
        rows.append(row)  # Add row to the list

    # Convert list of rows into a DataFrame
    df = pd.DataFrame(rows)

    # Ensure correct column order: filename, individuum, protocol, then species IDs
    df = df[['filename', 'individuum', 'protocol'] + sorted_species_ids]

    # Save the DataFrame as a CSV file
    output_path = os.path.join(results_folder, f"{species}{FILE_ENDING_KRAKEN_ALL_READS_COMBINED_ANALYSIS_CSV}")
    df.to_csv(output_path, index=False)
    print_info(f"Saved combined Kraken2 report to: {output_path}")

    # combine by individuum
    individuum_combined = df.groupby(['individuum']).sum(numeric_only=True).reset_index()
    individuum_combined_output_path = os.path.join(results_folder, f"{species}{FILE_ENDING_KRAKEN_BY_INDIVIDUAL_COMBINED_ANALYSIS_CSV}")
    individuum_combined.to_csv(individuum_combined_output_path, index=False)
    print_info(f"Saved combined Kraken2 report to: {individuum_combined_output_path}")

    # combine by protocol
    protocol_combined = df.groupby(['protocol']).sum(numeric_only=True).reset_index()
    protocol_combined_output_path = os.path.join(results_folder, f"{species}{FILE_ENDING_KRAKEN_BY_PROTOCOL_COMBINED_ANALYSIS_CSV}")
    protocol_combined.to_csv(protocol_combined_output_path, index=False)
    print_info(f"Saved combined Kraken2 report to: {protocol_combined_output_path}")
    

def run_Kraken_per_species(species: str):
    print_species_execution(f"Processing species: {species}")

    if not config.is_process_step_enabled(
        PipelineStages.RAW_READS_PROCESSING, 
        RawReadsProcessingSteps.CONTAMINATION_ANALYSIS,
        species,
        default=False
    ):
        print_skipping(f"Contamination analysis is not enabled for {species} in the config. Skipping this step.")
        return

    duplicates_removed_folder = get_folder_path_species_processed_duplicates_removed(species)
    
    fastq_files = get_files_in_folder_matching_pattern(
        duplicates_removed_folder,
        f"*{FILE_ENDING_DUPLICATES_REMOVED_FASTQ_GZ}"
    )

    # if no files are found, skip species
    if not fastq_files:
        print_warning(f"No duplicate removed FASTQ.GZ files found for species {species}.")
        return

    print_info(f"Found {len(fastq_files)} relevant FASTQ.GZ files for species {species}")
    print_debug(f"Files found: {fastq_files}")

    # Filter out files that contain "LB" or "EB" (Library Blanks or Extraction Blanks)
    fastq_files_filtered = [
        f for f in fastq_files
        if "LB" not in get_filename_from_path(f) and "EB" not in get_filename_from_path(f)
    ]

    excluded_files = [f for f in fastq_files if f not in fastq_files_filtered]
    if excluded_files:
        print_debug(f"Excluded {len(excluded_files)} files (LB/EB): {[get_filename_from_path(f) for f in excluded_files]}")

    # if no files remain after filtering, skip species
    if not fastq_files_filtered:
        print_warning(f"No FASTQ.GZ files found for species {species} after filtering.")
        return

    print_debug(f"Found {len(fastq_files_filtered)} relevant FASTQ.GZ files after filtering")
    print_debug(f"Files to process: {[get_filename_from_path(f) for f in fastq_files_filtered]}")

    # Run Kraken on each filtered file
    for fastq_file in fastq_files_filtered:

        filename_without_ext = get_filename_from_path_without_extension(fastq_file)
        processed_folder = get_folder_path_species_processed_qc_kraken(species)
        results_folder = get_folder_path_species_results_qc_kraken(species)

        # Define output file paths based on the input filename
        kraken_report_tsv = os.path.join(processed_folder, f"{filename_without_ext}{FILE_ENDING_KRAKEN_REPORT_TSV}")
        analysis_output_txt = os.path.join(results_folder, f"{filename_without_ext}{FILE_ENDING_KRAKEN_TOP5_ANALYSIS_TSV}") # Define analysis output path

        run_kraken_on_file(species, fastq_file, kraken_report_tsv)

        create_kraken_top5_analysis(kraken_report_tsv, analysis_output_txt)

    combine_kraken2_top5_analysis(species)

    print_info(f"Finished processing species: {species}")

def all_species_run_Kraken():
    print_step_execution("Starting Kraken analysis for all species on individual deduplicated FASTQ files.")

     # get the Kraken database path from the config
    Kraken_db = get_processing_settings(RawReadsProcessingSteps.CONTAMINATION_CHECK).get(ContaminationCheckSettings.KRAKEN_DB.value)

    # Ensure Kraken2 database path is set and exists
    if not Kraken_db:
        print_error("Kraken2 database path is not set. Please check your pipeline configuration.")
        return

    for species in FOLDER_SPECIES:
        run_Kraken_per_species(species)

    print_success("Finished Kraken analysis for all species.")