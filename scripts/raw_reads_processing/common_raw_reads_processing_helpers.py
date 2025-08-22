import os
import common_aDNA_scripts as common
from ref_genome_processing.common_ref_genome_processing_helpers import is_species_combined_reads_file_exists, is_species_individual_reads_file_exists, is_species_individual_and_combined_reads_file_exists, get_individual_from_file

def get_adapter_removed_path_for_paired_raw_reads(species, paired_read_file_path_list: list) -> str:
    filename_new = os.path.basename(paired_read_file_path_list[0]).replace("_R1_","_").replace("_R2_","_").replace(common.FILE_ENDING_FASTQ_GZ, common.FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ )
    return os.path.join(common.get_folder_path_species_processed_adapter_removed(species), filename_new)

def get_quality_filtered_path_for_adapter_removed_reads(species: str, adapter_removed_file_path: str) -> str:
    output_file = os.path.basename(adapter_removed_file_path).replace(common.FILE_ENDING_ADAPTER_REMOVED_FASTQ_GZ, common.FILE_ENDING_QUALITY_FILTERED_FASTQ_GZ)
    return os.path.join(common.get_folder_path_species_processed_quality_filtered(species), output_file)

def get_deduplication_path_for_quality_filtered_reads(species: str, quality_filtered_file_path: str) -> str:
    output_file = os.path.basename(quality_filtered_file_path).replace(common.FILE_ENDING_QUALITY_FILTERED_FASTQ_GZ, common.FILE_ENDING_DUPLICATES_REMOVED_FASTQ_GZ)
    return os.path.join(common.get_folder_path_species_processed_duplicates_removed(species), output_file)
