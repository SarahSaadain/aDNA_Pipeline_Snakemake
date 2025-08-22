import os
from common_aDNA_scripts import *
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

from ref_genome_processing.convert_mapped_sam2bam import execute_convert_sam_to_bam

def execute_bwa_map_mtDNA_to_refgenome(input_file_path:str, ref_genome_path:str, output_file_path:str, threads:int = THREADS_DEFAULT):
    
    print_info(f"Mapping {input_file_path} to reference genome ...")

    if not os.path.exists(input_file_path):
        raise Exception(f"Read file {input_file_path} does not exist!")
    
    if not os.path.exists(ref_genome_path):
        raise Exception(f"Reference genome file {ref_genome_path} does not exist!")
    
    if os.path.exists(output_file_path):
        print_info(f"Output file {output_file_path} already exists! Skipping!")
        return
    
    command_bwa = f"{PROGRAM_PATH_BWA} {PROGRAM_PATH_BWA_MEM} -M -T 50 -t {str(threads)} {ref_genome_path} {input_file_path} > {output_file_path}"
    print_debug(f"Executing command: {command_bwa}")

    try:
        subprocess.run(command_bwa, shell=True, check=True)
        print_success(f"Mapping {input_file_path} to reference genome complete")
    except Exception as e:
        print_error(f"Failed to run bwa for {input_file_path}: {e}")

def map_mtdna_to_refgenome_for_species(species: str):
    print_info(f"Mapping mtdna to reference genome for species {species} ...")

    #get reads
    read_folder = get_folder_path_species_raw_mtdna(species)
    list_of_mtrna_files = get_files_in_folder_matching_pattern(read_folder, f"*{FILE_ENDING_FASTA}")

    if len(list_of_mtrna_files) == 0:
        print_warning(f"No mtdna reads found for species {species}. Skipping.")
        return
    
    print_debug(f"Found {len(list_of_mtrna_files)} mtdna reads for species {species}")
    print_debug(f"mtdna reads: {list_of_mtrna_files}")

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_tuple in ref_genome_list:
        
        # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
        ref_genome_id = ref_genome_tuple[0]
        ref_genome_path = ref_genome_tuple[1]

        output_folder = get_folder_path_species_processed_refgenome_mtdna_mapped(species, ref_genome_id)
        

        for mtrna_read_file_path in list_of_mtrna_files:

            print_info(f"Mapping {mtrna_read_file_path} to reference genome {ref_genome_path} ...")

            sam_file_name = common_rgp.get_sam_file_name_for_read_file_and_ref_genome(mtrna_read_file_path, ref_genome_id)
            sam_file_path = os.path.join(output_folder, sam_file_name)

            bam_file_name = common_rgp.get_bam_file_name_for_sam_file(sam_file_path)
            sorted_bam_file_name = common_rgp.get_sorted_bam_file_name_for_bam_file(bam_file_name)

            bam_file_path = os.path.join(output_folder, bam_file_name)
            sorted_bam_file_path = os.path.join(output_folder, sorted_bam_file_name)

            try:

                if os.path.exists(sam_file_path):
                    print_info(f"Sam file {sam_file_path} already exists. Skipping mapping.")
                elif os.path.exists(sorted_bam_file_path):
                    print_info(f"Sorted BAM file {sorted_bam_file_path} already exists. Skipping mapping.")
                else: 
                    execute_bwa_map_mtDNA_to_refgenome(mtrna_read_file_path, ref_genome_path, sam_file_path)

                if not os.path.exists(sam_file_path):
                    print_error(f"Sam file {sam_file_path} was not created. Skipping conversion to BAM.")
                    continue

                print_info(f"Converting SAM file {sam_file_path} to BAM format ...")
                bam_file_name = common_rgp.get_bam_file_name_for_sam_file(sam_file_path)
                sorted_bam_file_name = common_rgp.get_sorted_bam_file_name_for_bam_file(bam_file_name)

                bam_file_path = os.path.join(output_folder, bam_file_name)
                sorted_bam_file_path = os.path.join(output_folder, sorted_bam_file_name)

                execute_convert_sam_to_bam(sam_file_path, bam_file_path, sorted_bam_file_path)
            except Exception as e:
                print_error(f"Failed to map {mtrna_read_file_path} to reference genome {ref_genome_path}: {e}")

    print_success(f"Mapping mtdna to reference genome for species {species} complete")

def all_species_map_mtdna_to_refgenome():
    print_step_execution("Mapping mtdna to reference genome for all species")

    for species in FOLDER_SPECIES: 
        map_mtdna_to_refgenome_for_species(species)

    print_info("Mapping mtdna to reference genome for all species complete")