import os
import subprocess

import common_aDNA_scripts as common_adna


def run_ecmsd_on_file(species: str, fastq_file_path: str, output_folder: str, threads: int = common_adna.THREADS_DEFAULT):
    common_adna.print_info(f"Running ECMSD on file: {common_adna.get_filename_from_path(fastq_file_path)}")

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Check for existing result file
    taxon_counts_output = os.path.join(output_folder, "Mito_summary.txt")
    if os.path.exists(taxon_counts_output):
        common_adna.print_skipping(f"Output file already exists: {taxon_counts_output}")
        return

    # Construct ECMSD command
    ecmsd_command = [
        "bash", common_adna.PROGRAM_PATH_ECMSD,
        "--fwd", fastq_file_path,
        "--out", output_folder,
        "--threads", str(threads),
        "--Binsize", "1000",
        "--RMUS-threshold", "0.15",
        "--mapping_quality", "20",
        "--taxonomic-hierarchy", "genus",
        "--force"
    ]

    common_adna.run_command(
        command=ecmsd_command,
        description=f"ECMSD on {common_adna.get_filename_from_path(fastq_file_path)}",
        cwd=common_adna.get_folder_path_resources_ecmsd()
    )
    common_adna.print_success(f"ECMSD completed for {common_adna.get_filename_from_path(fastq_file_path)}")

def run_ecmsd_per_species(species: str):
    common_adna.print_species_execution(f"Processing species: {species}")

    if not common_adna.config.is_process_step_enabled(
        common_adna.PipelineStages.RAW_READS_PROCESSING, 
        common_adna.RawReadsProcessingSteps.CONTAMINATION_ANALYSIS,
        species,
        default=False
    ):
        common_adna.print_skipping(f"Contamination analysis is not enabled for {species} in the config. Skipping this step.")
        return

   #get reads
    read_folder = common_adna.get_folder_path_species_processed_prepared_for_ref_genome(species)
    list_of_read_files = common_adna.get_files_in_folder_matching_pattern(read_folder, f"*{FILE_ENDING_FASTQ_GZ}")

    if len(list_of_read_files) == 0:
        common_adna.print_warning(f"No reads found for species {species}.")
        return
    
    common_adna.print_debug(f"Found {len(list_of_read_files)} read files for species {species}.")
    common_adna.print_debug(f"Read files: {list_of_read_files}")

    for fastq_file in list_of_read_files:
        filename_wo_ext = common_adna.get_filename_from_path_without_extension(fastq_file)
        output_folder = os.path.join(common_adna.get_folder_path_species_results_qc_ecmsd(species), filename_wo_ext)

        run_ecmsd_on_file(species, fastq_file, output_folder)


def all_species_run_ecmsd():
    common_adna.print_step_execution("Starting ECMSD analysis for all species on deduplicated FASTQ files.")

    for species in common_adna.FOLDER_SPECIES:
        run_ecmsd_per_species(species)

    common_adna.print_success("Finished ECMSD analysis for all species.")