import os
import csv

from Bio import SeqIO
from common_aDNA_scripts import *
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp


def check_extracted_region_for_species(species):

    print_info(f"Checking extracted region for species {species} ...")

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_tuple in ref_genome_list:

        print_info(f"Checking extracted region for reference genome {ref_genome_tuple[0]} for species {species}")
        
        # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
        ref_genome_id = ref_genome_tuple[0]
        #ref_genome_path = ref_genome_tuple[1]

        # Get the folder containing the extracted region sequences
        extracted_region_folder = get_folder_path_species_processed_refgenome_mtdna_extracted_sequence(species, ref_genome_id)

        # Get the list of extracted region files
        extracted_region_files = get_files_in_folder_matching_pattern(extracted_region_folder, f"*{FILE_ENDING_FASTA}")
        
        if len(extracted_region_files) == 0:
            print_warning(f"No extracted region files found for reference genome {ref_genome_id} for species {species}. Skipping.")
            continue
        
        print_debug(f"Found {len(extracted_region_files)} extracted region files for reference genome {ref_genome_id} for species {species}.")
        print_debug(f"Extracted region files: {extracted_region_files}")
        
        output_folder = get_folder_path_species_results_refgenome_mtdna_regions(species, ref_genome_id)
        
        output_file = os.path.join(output_folder, f"{species}_{ref_genome_id}_extracted_region_analysis.tsv")

        if os.path.exists(output_file):
            print_info(f"Output file {output_file} already exists. Skipping analysis.")
            continue
        
        results = []

        for file_path in extracted_region_files:

            print_info(f"Processing file {file_path} ...")

            filename = get_filename_from_path(file_path)

            # check if the file is empty
            if os.path.getsize(file_path) == 0:
                print_warning(f"File {file_path} is empty. Skipping.")
                continue

            with open(file_path, "r") as file:
                record = next(SeqIO.parse(file, "fasta"))  # Since each file has one read
                sequence = str(record.seq)
                
                total_length = len(sequence)
                non_n_count = sum(1 for base in sequence if base.upper() != "N")
                non_n_percentage = (non_n_count / total_length) * 100 if total_length > 0 else 0
                
                results.append([filename, total_length, non_n_count, round(non_n_percentage, 2)])

                # Print the results
                print_info(f"Results for {filename}:")
                print_info(f"Total length: {total_length}")
                print_info(f"Non-N count: {non_n_count}")
                print_info(f"Non-N percentage: {non_n_percentage:.2f}%")

        # Write results to a TSV file
        with open(output_file, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t")
            writer.writerow(["Filename", "Total Length", "Non-N Count", "Non-N Percentage"])
            writer.writerows(results)

    print_info(f"Consensus sequence of {species} created and mapped successfully.")

def all_species_check_extracted_region():
    print_step_execution("Check extracted regions for all species")

    for species in FOLDER_SPECIES: 
        check_extracted_region_for_species(species)

    print_info("All species extracted regions checked successfully.")