import os
import common_aDNA_scripts as common

def get_species_combined_read_path(species: str) -> str:
    output_folder = common.get_folder_path_species_processed_prepared_for_ref_genome(species)
    return os.path.join(output_folder, f"{species}_combined{common.FILE_ENDING_FASTQ_GZ}")

def is_species_combined_reads_file_exists(species: str) -> bool:
    combined_reads_path = get_species_combined_read_path(species)
    return os.path.exists(combined_reads_path)

def create_species_individual_combined_read_filepath(species: str, individual: str) -> str:
    output_folder = common.get_folder_path_species_processed_prepared_for_ref_genome(species)
    return os.path.join(output_folder, f"{individual}{common.FILE_ENDING_FASTQ_GZ}")

def is_species_individual_reads_file_exists(species: str, individual: str) -> bool:
    individual_reads_path = create_species_individual_combined_read_filepath(species, individual)
    return os.path.exists(individual_reads_path)

def is_species_individual_and_combined_reads_file_exists(species: str, individual: str) -> bool:
    individual_reads_exists = is_species_individual_reads_file_exists(species, individual)
    combined_reads_exists = is_species_combined_reads_file_exists(species)
    
    return individual_reads_exists and combined_reads_exists


def get_individual_from_file(file_path: str) -> str:

    filename = common.get_filename_from_path_without_extension(file_path)
    
    #the indivisual name is expected to be the first part of the filename, e.g. "individual_protocol_R1_001.fastq.gz"
    
    #check if filename has an underscore
    if '_' not in filename:
        raise ValueError(f"Filename '{filename}' does not contain an underscore to separate individual name.")

    parts = filename.split('_')

    return parts[0]  # Assuming the second part is the individual name

def get_sam_file_name_for_read_file_and_ref_genome(read_file_path: str, ref_genome_id: str) -> str:
    read_name = common.get_filename_from_path_without_extension(read_file_path)
    return f"{read_name}_{ref_genome_id}{common.FILE_ENDING_SAM}"

def get_sam_file_path_for_read_file_and_ref_genome(species: str, read_file_path: str, ref_genome_id: str) -> str:
    output_folder = common.get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)
    sam_file_name = get_sam_file_name_for_read_file_and_ref_genome(read_file_path, ref_genome_id)
    return os.path.join(output_folder, sam_file_name)

def get_bam_file_name_for_sam_file(sam_file_path: str) -> str:
    sam_file_name = common.get_filename_from_path_without_extension(sam_file_path)
    return f"{sam_file_name}{common.FILE_ENDING_BAM}"

def get_bam_file_path_for_sam_file(species: str, ref_genome_id: str, sam_file_path: str) -> str:
    output_folder = common.get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)
    bam_file_name = get_bam_file_name_for_sam_file(sam_file_path)
    return os.path.join(output_folder, bam_file_name)

def get_sorted_bam_file_name_for_bam_file(bam_file_path: str) -> str:
    bam_file_name = common.get_filename_from_path_without_extension(bam_file_path)
    return f"{bam_file_name}{common.FILE_ENDING_SORTED_BAM}"

def get_sorted_bam_file_path_for_bam_file(species: str, ref_genome_id: str, bam_file_path: str) -> str:
    output_folder = common.get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)
    sorted_bam_file_name = get_sorted_bam_file_name_for_bam_file(bam_file_path)
    return os.path.join(output_folder, sorted_bam_file_name)

def get_rescaled_bam_path_for_sorted_bam_path(bam_file_path: str) -> str:
    return bam_file_path.replace(common.FILE_ENDING_SORTED_BAM,common.FILE_ENDING_RESCALED_BAM)

def get_reference_genome_file_list_for_species(species: str) -> list[tuple[str, str]]:
 
    # get ref genome
    ref_genome_folder = common.get_folder_path_species_raw_ref_genome(species)

    # add fna files to reference genome list
    reference_genome_files = common.get_files_in_folder_matching_pattern(ref_genome_folder, f"*{common.FILE_ENDING_FNA}")
    # add fasta files to reference genome list
    reference_genome_files += common.get_files_in_folder_matching_pattern(ref_genome_folder, f"*{common.FILE_ENDING_FASTA}")

    # add fa files to reference genome list
    reference_genome_files += common.get_files_in_folder_matching_pattern(ref_genome_folder, f"*{common.FILE_ENDING_FA}")

    if len(reference_genome_files) == 0:
        raise Exception(f"No reference genome found for species {species}.")
    
    common.print_debug(f"Found {len(reference_genome_files)} reference genome files for species {species}.")
    common.print_debug(f"Reference genome files: {reference_genome_files}")

    # return as tuple of (filename without extension, filepath)
    reference_genome_files_with_filename = [(os.path.splitext(os.path.basename(f))[0], f) for f in reference_genome_files]

    return reference_genome_files_with_filename