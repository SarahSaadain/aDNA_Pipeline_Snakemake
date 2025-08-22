import common_aDNA_scripts as common_adna

import additional_analysis.mtdna_analysis.determine_mtdna_step1_map_to_ref_genome as determine_mtdna_step1_map_to_ref_genome
import additional_analysis.mtdna_analysis.determine_mtdna_step2_determine_regions as determine_mtdna_step2_determine_regions
import additional_analysis.mtdna_analysis.determine_mtdna_step4_extract_coi_regions as determine_mtdna_step4_extract_coi_regions
import additional_analysis.mtdna_analysis.determine_mtdna_step3_create_and_map_consensus_sequence as determine_mtdna_step3_create_and_map_consensus_sequence
import additional_analysis.mtdna_analysis.determine_mtdna_step5_check_extracted_regions_for_content as determine_mtdna_step5_check_extracted_regions_for_content

def pipeline_mtdna_analysis():

    common_adna.print_step_execution("Starting mtDNA analysis pipeline ...")

    if not common_adna.config.is_process_step_enabled(
        common_adna.PipelineStages.POST_PROCESSING, 
        common_adna.PostProcessingSteps.MTDNA_ANALYSIS.value
    ):
        common_adna.print_skipping("mapDamage analysis is not enabled in the config. Skipping this step.")
        return

    determine_mtdna_step1_map_to_ref_genome.all_species_map_mtdna_to_refgenome()
    determine_mtdna_step2_determine_regions.all_species_mtdna_get_regions()
    determine_mtdna_step3_create_and_map_consensus_sequence.all_species_create_consensus_sequence()
    determine_mtdna_step4_extract_coi_regions.all_species_extract_mtdna_region()
    determine_mtdna_step5_check_extracted_regions_for_content.all_species_check_extracted_region()