import os
from common_aDNA_scripts import *
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp


def execute_samtools_get_read_regions(bam_file: str, output_file: str, threads: int=THREADS_DEFAULT):

    if not os.path.exists(bam_file):
        raise Exception(f"BAM file {bam_file} does not exist.")
    
    if os.path.exists(output_file):
        print_warning(f"Output file {output_file} already exists. Skipping.")
        return
    
    # Exclude secondary alignments (SAM flag 0x100) to ensure only primary alignments 
    # are passed to bedtools for BED conversion
    command = (
        f"{PROGRAM_PATH_SAMTOOLS} {PROGRAM_PATH_SAMTOOLS_VIEW} -h -@ {threads} -F 0x100 {bam_file} | "
        f"{PROGRAM_PATH_BEDTOOLS} {PROGRAM_PATH_BAMTOBED} -i > {output_file}"
    )
    print_debug(f"Executing command: {command}")
    
    try:
      # Execute the command
        subprocess.run(command, shell=True, check=True)
        print_success(f"Regions for {bam_file} have been written to {output_file}")
    except Exception as e:
        print_error(f"Failed to extract regions for {bam_file}: {e}")

def mtdna_get_regions_for_species(species):
    print_info(f"Determining mtdna regions for species: {species}")

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_tuple in ref_genome_list:
        
        # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
        ref_genome_id = ref_genome_tuple[0]
        #ref_genome_path = ref_genome_tuple[1]
    
        mapped_folder = get_folder_path_species_processed_refgenome_mtdna_mapped(species, ref_genome_id)
        
        bam_files = get_files_in_folder_matching_pattern(mapped_folder, f"*{FILE_ENDING_SORTED_BAM}")

        if len(bam_files) == 0:
            print_warning(f"No BAM files found for species {species}. Skipping.")
            continue
        
        print_debug(f"Found {len(bam_files)} BAM files for species {species}.")
        print_debug(f"BAM files: {bam_files}")
        
        for bam_file in bam_files:
            print_info(f"Determining mtdna regions for {bam_file}")

            bam_file_name_wo_ext = get_filename_from_path_without_extension(bam_file)

            result_folder = get_folder_path_species_results_refgenome_mtdna_regions(species, ref_genome_id)
            result_file_path = os.path.join(result_folder, f"{bam_file_name_wo_ext}_mtdna_region{FILE_ENDING_BED}")

            if os.path.exists(result_file_path):
                print_info(f"Result file {result_file_path} already exists for species {species}. Skipping.")
                continue
            
            execute_samtools_get_read_regions(bam_file, result_file_path)

    print_info(f"Finished determining mtdna regions for species {species}")
    
   
def all_species_mtdna_get_regions():

    print_step_execution("Determining mtdna regions for all species")
    for species in FOLDER_SPECIES: 
        mtdna_get_regions_for_species(species)
    print_info("Finished determining mtdna regions for all species")