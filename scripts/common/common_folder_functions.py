import os
from common.common_constants import *
from common.common_config import *

#####################
# Folder paths
#####################

def is_species_folder(folder_name: str) -> bool:
    return folder_name in [config_dict['species'][s]['folder_name'] for s in config_dict['species']]

def check_folder_exists_or_create(folder_path: str) -> None:
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

def get_folder_aDNA() -> str:
    return PATH_ADNA_PROJECT

def get_folder_path_results() -> str:
    path = os.path.join(get_folder_aDNA(), FOLDER_RESULTS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_logs() -> str:
    path = os.path.join(get_folder_aDNA(), FOLDER_LOGS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_resources() -> str:
    path = os.path.join(get_folder_path_resources(), FOLDER_RESOURCES)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_resources_ecmsd() -> str:
    path = os.path.join(get_folder_path_resources(), FOLDER_ECMSD)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_results_plots() -> str:
    path = os.path.join(get_folder_path_results(), FOLDER_PLOTS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_scripts() -> str:
    path = os.path.join(get_folder_aDNA(), FOLDER_SCRIPTS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_scripts_plots(processing_folder: str) -> str:
    path = os.path.join(get_folder_path_scripts(), processing_folder, FOLDER_ANALYSIS, FOLDER_PLOTS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_raw(species: str) -> str:
    if not is_species_folder(species):
        raise Exception(f"Invalid species folder: {species}")

    path = os.path.join(get_folder_aDNA(), species, FOLDER_RAW)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species(species: str) -> str:
    if not is_species_folder(species):
        raise Exception(f"Invalid species folder: {species}")

    path = os.path.join(get_folder_aDNA(), species)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed(species: str) -> str:
    path = os.path.join(get_folder_path_species(species), FOLDER_PROCESSED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results(species: str) -> str:
    path = os.path.join(get_folder_path_species(species), FOLDER_RESULTS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_scripts(species: str) -> str:
    path = os.path.join(get_folder_path_species(species), FOLDER_SCRIPTS)
    check_folder_exists_or_create(path)
    return path


def get_folder_path_species_logs(species: str) -> str:
    path = os.path.join(get_folder_path_species(species), FOLDER_LOGS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_resources(species: str) -> str:
    path = os.path.join(get_folder_path_species(species), FOLDER_RESOURCES)
    check_folder_exists_or_create(path)
    return path



def get_folder_path_species_raw_reads(species: str) -> str:
    path = os.path.join(get_folder_path_species_raw(species), FOLDER_READS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_raw_mtdna(species: str) -> str:
    path = os.path.join(get_folder_path_species_raw(species), FOLDER_MTDNA)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_raw_ref_genome(species: str) -> str:
    path = os.path.join(get_folder_path_species_raw(species), FOLDER_REFERENCE_GENOMES)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), reference_genome_id)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_mapped(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome(species, reference_genome_id), FOLDER_MAPPED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_consensus_sequences(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome(species, reference_genome_id), FOLDER_CONSENSUS_SEQUENCES)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_mtdna(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome(species, reference_genome_id), FOLDER_MTDNA)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_mtdna_extracted_sequence(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome_mtdna(species, reference_genome_id), FOLDER_EXTRACTED_SEQUENCES)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_mtdna_mapped(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome_mtdna(species, reference_genome_id), FOLDER_MAPPED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_mtdna_consensus_sequences(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome_mtdna(species, reference_genome_id), FOLDER_CONSENSUS_SEQUENCES)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_mtdna_consensus_sequences_mapped(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome_mtdna(species, reference_genome_id), FOLDER_CONSENSUS_SEQUENCES_MAPPED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_adapter_removed(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_ADAPTER_REMOVED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_quality_filtered(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_QUALITY_FILTERED)
    check_folder_exists_or_create(path)
    return path 

def get_folder_path_species_processed_duplicates_removed(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_DUPLICATES_REMOVED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_qc(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_QUALITYCONTROL)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_qc_reads_processing(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed_qc(species), FOLDER_PROCESSED_READS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_qc_read_length_distribution(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed_qc(species), FOLDER_READ_LENGTH_DISTRIBUTION)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_qc_centrifuge(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed_qc(species), FOLDER_CENTRIFUGE)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_qc_kraken(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed_qc(species), FOLDER_KRAKEN)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_prepared_for_ref_genome(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_PREPARED_FOR_REF_GENOME)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_genomedelta(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_GENOMEDELTA)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_endogenous_reads(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome(species, reference_genome_id), FOLDER_ENDOGENOUS_READS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_concatenated(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_CONCATENATED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_non_concatenated(species: str) -> str:
    path = os.path.join(get_folder_path_species_processed(species), FOLDER_NON_CONCATENATED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_coverage(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome(species, reference_genome_id), FOLDER_DEPTH_BREADTH)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_damage(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome(species, reference_genome_id), FOLDER_DAMAGE_ANALYSIS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_processed_refgenome_damage_individual(species: str, reference_genome_id: str, individual_id: str) -> str:
    path = os.path.join(get_folder_path_species_processed_refgenome_damage(species, reference_genome_id), individual_id)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_genomedelta(species: str) -> str:
    path = os.path.join(get_folder_path_species_results(species), FOLDER_GENOMEDELTA)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results(species), reference_genome_id)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_mtdna(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome(species, reference_genome_id), FOLDER_MTDNA)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_damage(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome(species, reference_genome_id), FOLDER_DAMAGE_ANALYSIS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_damage_individual(species: str, reference_genome_id: str, individual_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_damage(species, reference_genome_id), individual_id)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_mtdna_regions(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_mtdna(species, reference_genome_id), FOLDER_REGIONS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc(species: str) -> str:
    path = os.path.join(get_folder_path_species_results(species), FOLDER_QUALITYCONTROL)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_centrifuge(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc(species), FOLDER_CENTRIFUGE)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_ecmsd(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc(species), FOLDER_ECMSD)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_kraken(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc(species), FOLDER_KRAKEN)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_fastqc(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc(species), FOLDER_FASTQC)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_fastqc_raw(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_fastqc(species), FOLDER_RAW)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_fastqc_adapter_removed(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_fastqc(species), FOLDER_ADAPTER_REMOVED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_fastqc_quality_filtered(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_fastqc(species), FOLDER_QUALITY_FILTERED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_fastqc_duplicates_removed(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_fastqc(species), FOLDER_DUPLICATES_REMOVED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_multiqc(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc(species), FOLDER_MULTIQC)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_multiqc_raw(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_multiqc(species), FOLDER_RAW)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_multiqc_adapter_removed(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_multiqc(species), FOLDER_ADAPTER_REMOVED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_multiqc_quality_filtered(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_multiqc(species), FOLDER_QUALITY_FILTERED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_multiqc_duplicates_removed(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc_multiqc(species), FOLDER_DUPLICATES_REMOVED)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_coverage(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome(species, reference_genome_id), FOLDER_DEPTH_BREADTH)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_mitochondria(species: str) -> str:
    path = os.path.join(get_folder_path_species_results(species), FOLDER_MITOCHONDRIA)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_plots(species: str) -> str:
    path = os.path.join(get_folder_path_species_results(species), FOLDER_PLOTS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_plots(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome(species, reference_genome_id), FOLDER_PLOTS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_plots_reads_processing(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_plots(species), FOLDER_PROCESSED_READS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_plots_depth(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_plots(species, reference_genome_id), FOLDER_DEPTH)
    check_folder_exists_or_create(path) 
    return path

def get_folder_path_species_results_refgenome_plots_depth_sample(species: str, reference_genome_id: str, sample_name: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_plots_depth(species, reference_genome_id), sample_name)
    check_folder_exists_or_create(path) 
    return path

def get_folder_path_species_results_refgenome_plots_breadth(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_plots(species, reference_genome_id), FOLDER_BREADTH)
    check_folder_exists_or_create(path) 
    return path

def get_folder_path_species_results_refgenome_plots_breadth_sample(species: str, reference_genome_id: str, sample_name: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_plots_breadth(species, reference_genome_id), sample_name)
    check_folder_exists_or_create(path) 
    return path

def get_folder_path_species_results_refgenome_plots_endogenous_reads(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_plots(species, reference_genome_id), FOLDER_ENDOGENOUS_READS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_plots_snp(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome_plots(species, reference_genome_id), FOLDER_SNP)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_plots_read_length_distribution(species: str) -> str:    
    path = os.path.join(get_folder_path_species_results_plots(species), FOLDER_READ_LENGTH_DISTRIBUTION)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_plots_contamination(species: str) -> str:    
    path = os.path.join(get_folder_path_species_results_plots(species), FOLDER_CONTAMINATION)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_refgenome_special_sequences(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome(species, reference_genome_id), FOLDER_SPECIAL_SEQUENCES, reference_genome_id)
    check_folder_exists_or_create(path) 
    return path

def get_folder_path_species_results_refgenome_endogenous_reads(species: str, reference_genome_id: str) -> str:
    path = os.path.join(get_folder_path_species_results_refgenome(species, reference_genome_id), FOLDER_ENDOGENOUS_READS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_reads_processing(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc(species), FOLDER_PROCESSED_READS)
    check_folder_exists_or_create(path)
    return path

def get_folder_path_species_results_qc_read_length_distribution(species: str) -> str:
    path = os.path.join(get_folder_path_species_results_qc(species), FOLDER_READ_LENGTH_DISTRIBUTION)
    check_folder_exists_or_create(path)
    return path

