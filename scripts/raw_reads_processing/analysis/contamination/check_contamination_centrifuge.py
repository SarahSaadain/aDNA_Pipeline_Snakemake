import os
import subprocess

from common_aDNA_scripts import *


def run_centrifuge_on_file(species: str, fastq_file_path: str, centrifuge_output_txt: str, centrifuge_report_tsv: str, threads: int = THREADS_DEFAULT):
   
    print_info(f"Running Centrifuge on file: {get_filename_from_path(fastq_file_path)}")

             # Check if output files already exist
    if os.path.exists(centrifuge_output_txt) and os.path.exists(centrifuge_report_tsv):
        print_skipping(f"Output files for {get_filename_from_path(fastq_file_path)} already exist.")
        return

    # get the Centrifuge database path from the config
    centrifuge_db = config.get_pipeline_settings(RawReadsProcessingSteps.CONTAMINATION_CHECK).get(ContaminationCheckSettings.CENTRIFUGE_DB.value)

    if not centrifuge_db:
        print_error("Centrifuge database path is not set. Please check your configuration.")
        return

    # https://ccb.jhu.edu/software/centrifuge/manual.shtml#usage

    # Construct the centrifuge command
    centrifuge_command = [
        PROGRAM_PATH_CENTRIFUGE,
        "-x", centrifuge_db,
        "-U", fastq_file_path,
        "-S", centrifuge_output_txt,
        "--report-file", centrifuge_report_tsv,
        "--threads", str(threads), # Number of threads        
        "--seed", "999"
    ]

    print_debug(f"Centrifuge command: {' '.join(centrifuge_command)}")

    # Execute the command
    try:
        result = subprocess.run(centrifuge_command, check=True, capture_output=True, text=True)
        print_success(f"Centrifuge analysis complete for {get_filename_from_path(fastq_file_path)}")
        # Optionally print stdout/stderr for debugging
        print_debug("Centrifuge stdout:\n" + result.stdout)
        print_debug("Centrifuge stderr:\n" + result.stderr)
    except subprocess.CalledProcessError as e:
        print_error(f"Centrifuge failed for {get_filename_from_path(fastq_file_path)} with error: {e}")
        print_error("Centrifuge stdout:\n" + e.stdout)
        print_error("Centrifuge stderr:\n" + e.stderr)
    except FileNotFoundError:
         print_error("Centrifuge command not found. Make sure Centrifuge is installed and in your PATH.")
    except Exception as e:
        print_error(f"An unexpected error occurred while running Centrifuge on {get_filename_from_path(fastq_file_path)}: {e}")

def analyze_centrifuge_output(output_file_path: str, taxon_counts_output_path: str):

    print_info(f"Analyzing Centrifuge output for taxon counts: {get_filename_from_path(output_file_path)}")

    # Check if output file already exists
    if os.path.exists(taxon_counts_output_path):
        print_skipping(f"Taxon counts output file {get_filename_from_path(taxon_counts_output_path)} already exists.")
        return

    # Construct the shell command pipeline
    # awk '$3 != 0 {print $3}': Filters lines where the 3rd column (abundance) is not 0 and prints the 3rd column (taxon ID)
    # sort: Sorts the taxon IDs
    # uniq -c: Counts occurrences of each unique taxon ID
    # sort -nr: Sorts the counts in reverse numerical order
    analysis_command = f"awk '$3 != 0 {{print $3}}' {output_file_path} | sort | uniq -c | sort -nr > {taxon_counts_output_path}"

    print_debug(f"Analysis command: {analysis_command}")

    # Execute the command
    try:
        # Use shell=True because we are using a pipeline with pipes (|) and redirection (>)
        subprocess.run(analysis_command, shell=True, check=True, capture_output=True, text=True)
        print_success(f"Taxon counts analysis complete. Results written to {get_filename_from_path(taxon_counts_output_path)}")
        # Optionally print stdout/stderr for debugging
        # print_debug("Analysis stdout:\n" + result.stdout)
        # print_debug("Analysis stderr:\n" + result.stderr)
    except subprocess.CalledProcessError as e:
        print_error(f"Taxon counts analysis failed for {get_filename_from_path(output_file_path)} with error: {e.returncode}")
        print_error("Analysis stdout:\n" + e.stdout)
        print_error("Analysis stderr:\n" + e.stderr)
    except Exception as e:
        print_error(f"An unexpected error occurred during taxon counts analysis of {get_filename_from_path(output_file_path)}: {e}")


def run_centrifuge_per_species(species: str):
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
        print_debug(f"Excluded {len(excluded_files)} files (LB/EB): {[os.path.basename(f) for f in excluded_files]}")

    # if no files remain after filtering, skip species
    if not fastq_files_filtered:
        print_warning(f"No FASTQ.GZ files found for species {species} after filtering.")
        return

    print_debug(f"Found {len(fastq_files_filtered)} relevant FASTQ.GZ files after filtering")
    print_debug(f"Files to process: {[os.path.basename(f) for f in fastq_files_filtered]}")

    # Run centrifuge on each filtered file
    for fastq_file in fastq_files_filtered:

        filename_without_ext = get_filename_from_path_without_extension(fastq_file)
        processed_folder = get_folder_path_species_processed_qc_centrifuge(species)
        results_folder = get_folder_path_species_results_qc_centrifuge(species)

        # Define output file paths based on the input filename
        centrifuge_output_txt = os.path.join(processed_folder, f"{filename_without_ext}{FILE_ENDING_CENTRIFUGE_OUTPUT_TXT}")
        centrifuge_report_tsv = os.path.join(processed_folder, f"{filename_without_ext}{FILE_ENDING_CENTRIFUGE_REPORT_TSV}")

        taxon_counts_output_txt = os.path.join(results_folder, f"{filename_without_ext}{FILE_ENDING_CENTRIFUGE_TAXON_COUNTS_TXT}") # Define analysis output path

        run_centrifuge_on_file(species, fastq_file, centrifuge_output_txt , centrifuge_report_tsv)

        analyze_centrifuge_output(centrifuge_output_txt, taxon_counts_output_txt)

    print_info(f"Finished processing species: {species}")

def all_species_run_centrifuge():
    print_step_execution("Starting Centrifuge analysis for all species on individual deduplicated FASTQ files.")

    # get the Centrifuge database path from the config
    centrifuge_db = get_processing_settings(RawReadsProcessingSteps.CONTAMINATION_CHECK).get(ContaminationCheckSettings.CENTRIFUGE_DB.value)

    if not centrifuge_db:
        print_error("Centrifuge database path is not set. Please check your configuration.")
        return

    for species in FOLDER_SPECIES:
        run_centrifuge_per_species(species)

    print_success("Finished Centrifuge analysis for all species.")