import os
from common_aDNA_scripts import *

import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

def all_species_prepare_ref_genome():

    print_step_execution("Preparing reference genomes for mapping for all species")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.PREPARE_REFERENCE_GENOME
    ):
        print_skipping("Preparing reference genomes for mapping is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES:
        species_prepare_ref_genome(species)
    
    print_success("Finished preparing reference genomes for mapping for all species")

def species_prepare_ref_genome(species: str):
    print_species_execution(f"Preparing reference genome(s) for {species} for mapping")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.PREPARE_REFERENCE_GENOME,
        species
    ):
        print_skipping(f"Preparing reference genomes for mapping is not enabled in the config for species {species}. Skipping this step.")
        return

    try:
        ref_genome_files = common_rgp.get_reference_genome_file_list_for_species(species)

        number_of_ref_genome_files = len(ref_genome_files)
        count_current = 0

        for ref_genome in ref_genome_files:

            count_current += 1

            print_info(f"[{count_current}/{number_of_ref_genome_files}] Indexing reference genome {ref_genome} ...")

            # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
            ref_genome_file_path = ref_genome[1]
            index_fasta_file_bwa(ref_genome_file_path)
            index_fasta_file_samtools(ref_genome_file_path)

        print_info(f"Finished preparing reference genome for {species} for mapping")

    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return
