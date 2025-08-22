import os
from common_aDNA_scripts import *

import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

def execute_samtools_count_mapped_reads(bam_file: str, threads: int=THREADS_DEFAULT) -> int:
    """Counts the number of mapped reads in a BAM file."""
    try:
        command = [PROGRAM_PATH_SAMTOOLS, PROGRAM_PATH_SAMTOOLS_VIEW, "-@", str(threads), "-c", "-F", "4", bam_file]

        result = run_command(
            command, 
            description=f"Counting mapped reads in {bam_file}",
            throw_error=True
        )

        return int(result)
    except Exception as e:
        print_error(f"Failed to count mapped reads in {bam_file}: {e}")
    return 0
def execute_samtools_count_total_reads(bam_file: str, threads: int=THREADS_DEFAULT) -> int:
    """Counts the total number of reads in a BAM file (mapped and unmapped)."""
    try:
        command = [PROGRAM_PATH_SAMTOOLS, PROGRAM_PATH_SAMTOOLS_VIEW, "-@", str(threads), "-c", bam_file]
        
        result = run_command(
            command, 
            description=f"Counting total reads in {bam_file}",
            throw_error=True
        )
       
        return int(result)
    except Exception as e:
        print_error(f"Failed to count total reads in {bam_file}: {e}")
    return 0


def determine_endogenous_reads_for_bam_file(bam_file: str, processed_folder: str):
    """Calculates endogenous reads for a single BAM file and saves to a separate file."""
    bam_filename = get_filename_from_path_without_extension(bam_file)

    output_filename = f"{bam_filename}{FILE_ENDING_ENDOGENOUS_READS_CSV}"
    target_file_path = os.path.join(processed_folder, output_filename)

    if os.path.exists(target_file_path):
        print_skipping(f"Result file already exists for {bam_filename}.")
        return

    print_info(f"Determining endogenous reads for {bam_filename} ...")

    proportion = 0.0

    print_debug(f"Counting mapped reads in {bam_filename} ...")
    mapped_reads = execute_samtools_count_mapped_reads(bam_file)

    print_debug(f"Counting total reads in {bam_filename} ...")
    total_reads = execute_samtools_count_total_reads(bam_file)

    if total_reads != 0:
        proportion = mapped_reads / total_reads

    print_info(f"Endogenous reads for {bam_filename}: {mapped_reads}/{total_reads} -> {proportion:.4f}")

    print_debug(f"Writing results for {bam_filename} to {target_file_path}")
    try:
        with open(target_file_path, "w") as result_file:
            
            result_file.write("Filename,MappedReads,TotalReads,Proportion\n")
            result_file.write(f"{bam_filename},{mapped_reads},{total_reads},{proportion}\n")
        
        print_debug(f"Wrote results for {bam_filename} to {target_file_path}")
    except IOError as e:
        print_error(f"Failed to write result file for {bam_filename} at {target_file_path}: {e}")


def combine_endogenous_reads_files(species: str, ref_genome_id: str):
    """Combines individual endogenous read files for a species into a single summary file."""
    
    print_info(f"Combining endogenous reads files for species: {species}")
    
    processed_folder = get_folder_path_species_processed_refgenome_endogenous_reads(species, ref_genome_id)
    result_folder = get_folder_path_species_results_refgenome_endogenous_reads(species, ref_genome_id)
    
    combined_file_path = os.path.join(result_folder, f"{species}{FILE_ENDING_ENDOGENOUS_READS_CSV}")

    # Check if there are any individual files to combine
    individual_files = get_files_in_folder_matching_pattern(processed_folder, f"*{FILE_ENDING_ENDOGENOUS_READS_CSV}")

    if not individual_files:
        print_warning(f"No individual endogenous read files found to combine for species {species}.")
        return

    print_debug(f"Found {len(individual_files)} individual endogenous read files to combine.")
    print_debug(f"Individual endogenous read files: {individual_files}")

    try:
        print_debug(f"Writing combined results to {combined_file_path}")
        with open(combined_file_path, "w") as combined_file:
            # Write header for the combined file
            combined_file.write("Filename,MappedReads,TotalReads,Proportion\n")

            # Sort files alphabetically for consistent output order
            individual_files.sort()

            for individual_file in individual_files:
                try:

                    print_debug(f"Reading individual file {individual_file}")
                    with open(individual_file, "r") as f_in:
                        # The individual file has a header. Skip it
                        header = f_in.readline()

                        # Write the rest of the file
                        data_line = f_in.readline().strip()
                        
                        # Simple check to avoid writing empty lines if files were corrupted/empty
                        if data_line:
                            combined_file.write(data_line + "\n")
                except IOError as e:
                     print_error(f"Could not read individual file {individual_file}: {e}")
                except Exception as e:
                    print_error(f"An unexpected error occurred reading {individual_file}: {e}")

        print_success(f"Successfully created combined endogenous reads file: {combined_file_path}")

    except IOError as e:
        print_error(f"Failed to write combined result file at {combined_file_path}: {e}")
    except Exception as e:
        print_error(f"An unexpected error occurred during combining files for {species}: {e}")


def determine_endogenous_reads_for_species(species: str, ref_genome_id: str):
    """Processes all BAM files for a species, creating individual result files."""
    print_info(f"Processing BAM files to determine endogenous reads for species: {species}")

    mapped_folder = get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)
    processed_folder = get_folder_path_species_processed_refgenome_endogenous_reads(species, ref_genome_id)

    bam_files = get_files_in_folder_matching_pattern(mapped_folder, f"*{FILE_ENDING_SORTED_BAM}")

    if len(bam_files) == 0:
        print_warning(f"No BAM files found in {mapped_folder} for species {species}.")
        return

    print_debug(f"Found {len(bam_files)} BAM files for species {species}: {bam_files}")

    for bam_file in bam_files:
        determine_endogenous_reads_for_bam_file(bam_file, processed_folder)

    print_info(f"Finished processing individual BAM files for species {species}.")

def all_species_determine_endogenous_reads():

    print_step_execution("Determining endogenous reads for all species")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.ENDOGENOUS_READS_ANALYSIS
    ):
        print_skipping("Determining endogenous reads is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 

        print_species_execution(f"Determining endogenous reads for species {species}")

        if not config.is_process_step_enabled(
            PipelineStages.REFERENCE_GENOME_PROCESSING, 
            ReferenceGenomeProcessingSteps.ENDOGENOUS_READS_ANALYSIS,
            species
        ):
            print_skipping(f"Creating consensus sequence for species {species} is disabled in the config. Skipping this species.")
            return

        try:
            ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
        except Exception as e:
            print_error(f"Failed to get reference genome files for species {species}: {e}")
            continue

        for ref_genome_tuple in ref_genome_list:

            # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
            ref_genome_id = ref_genome_tuple[0]
            #ref_genome_path = ref_genome_tuple[1]

            print_info(f"Determining endogenous reads for species {species} and reference genome {ref_genome_id}")

            # First, process each individual BAM file for the species
            determine_endogenous_reads_for_species(species, ref_genome_id)
            # Then, combine the individual results for this species
            combine_endogenous_reads_files(species, ref_genome_id)

    print_success("Endogenous reads determination completed for all species.")
