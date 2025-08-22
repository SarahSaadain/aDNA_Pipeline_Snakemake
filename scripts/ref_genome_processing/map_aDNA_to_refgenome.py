import os
import subprocess
from common_aDNA_scripts import *

import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp
import ref_genome_processing.convert_mapped_sam2bam as convert_sam2bam

def execute_bwa_map_aDNA_to_refgenome(input_file_path:str, ref_genome_path:str, output_file_path:str, threads:int = THREADS_DEFAULT):
    
    print_info(f"Mapping {input_file_path} to reference genome ...")

    if not os.path.exists(input_file_path):
        raise Exception(f"Read file {input_file_path} does not exist!")
    
    if not os.path.exists(ref_genome_path):
        raise Exception(f"Reference genome file {ref_genome_path} does not exist!")
    
    if os.path.exists(output_file_path):
        print_skipping(f"Output file {output_file_path} already exists!")
        return
    
    command_bwa = f"{PROGRAM_PATH_BWA} {PROGRAM_PATH_BWA_MEM} -t {str(threads)} {ref_genome_path} {input_file_path} > {output_file_path}"
    print_debug(f"BWA command: {command_bwa}")

    try:
        subprocess.run(command_bwa, shell=True, check=True)
        print_success(f"Mapping {input_file_path} to reference genome complete")
    except Exception as e:
        print_error(f"Failed to run bwa for {input_file_path}: {e}")

def map_aDNA_to_refgenome_for_species(species: str):
    print_species_execution(f"Mapping aDNA to reference genome for species {species} ...")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.MAP_READS_TO_REFERENCE_GENOME,
        species
    ):
        print_skipping(f"Mapping aDNA to reference genome is not enabled in the config for species {species}. Skipping this step.")
        return

    #get reads
    read_folder = get_folder_path_species_processed_prepared_for_ref_genome(species)
    list_of_read_files = get_files_in_folder_matching_pattern(read_folder, f"*{FILE_ENDING_FASTQ_GZ}")

    if len(list_of_read_files) == 0:
        print_warning(f"No reads found for species {species}.")
        return
    
    print_debug(f"Found {len(list_of_read_files)} read files for species {species}.")
    print_debug(f"Read files: {list_of_read_files}")

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return
    
    print_debug(f"Found {len(ref_genome_list)} reference genome files for species {species}.")
    print_debug(f"Reference genome files: {ref_genome_list}")
    
    number_of_entries = len(ref_genome_list) * len(list_of_read_files)
    count_current = 0

    for ref_genome_tuple in ref_genome_list:

        # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
        ref_genome_id = ref_genome_tuple[0]
        ref_genome_path = ref_genome_tuple[1]

        print_debug(f"Reference genome file: {ref_genome_id}")

        for read_file_path in list_of_read_files:

            count_current += 1

            print_info(f"[{count_current}/{number_of_entries}] Mapping {read_file_path} to reference genome {ref_genome_path} ...")

            sam_file_path = common_rgp.get_sam_file_path_for_read_file_and_ref_genome(species, read_file_path, ref_genome_id)

            print_debug(f"Output file: {sam_file_path}")
            
            # we only need to map if the sorted bam file does not exist
            bam_file_path = common_rgp.get_bam_file_path_for_sam_file(species, ref_genome_id, sam_file_path)
            sorted_bam_file_path = common_rgp.get_sorted_bam_file_path_for_bam_file(species, ref_genome_id, bam_file_path)

            print_debug(f"Sorted BAM file path: {sorted_bam_file_path}")

            if os.path.exists(sorted_bam_file_path):
                print_skipping(f"Sorted BAM file {sorted_bam_file_path} already exists!")
                continue

            if os.path.exists(sam_file_path):
                print_skipping(f"SAM file {sam_file_path} already exists!")
            else:
                execute_bwa_map_aDNA_to_refgenome(read_file_path, ref_genome_path, sam_file_path, THREADS_DEFAULT)

            # Convert SAM to BAM and sort
            # Add this here so the SAM to BAM conversion is done directly after mapping
            # This will help to reduce the space used by the SAM files as they can be very large and are not needed after conversion
            # if this step will be called later, we will require more space as first all SAM files will be created and then converted to BAM files
            convert_sam2bam.execute_convert_sam_to_bam(sam_file_path, bam_file_path, sorted_bam_file_path)


    print_success(f"Mapping aDNA to reference genome for species {species} complete")

def all_species_map_aDNA_to_refgenome():

    print_step_execution("Mapping aDNA to reference genome for all species")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.MAP_READS_TO_REFERENCE_GENOME
    ):
        print_skipping("Mapping aDNA to reference genome is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        map_aDNA_to_refgenome_for_species(species)

    print_info("Mapping aDNA to reference genome for all species complete")
