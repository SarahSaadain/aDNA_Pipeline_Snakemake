import pysam
import os
from collections import Counter
from common_aDNA_scripts import *
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

DEPTH_THRESHOLD = 1000
MINIMUM_SEQUENCE_LENGTH = 1000
MAXIMUM_SEQUENCE_LENGTH = 5000

def write_fasta_entry(special_reads_file_content: str, scaffold: str, start: int, end: int, sequence: str, depth_values: list, minimum_sequence_length: int=MINIMUM_SEQUENCE_LENGTH):
    """Writes a sequence entry to the output FASTA file."""
    if sequence:

        if len(sequence) < minimum_sequence_length:
            return

        avg_depth = round(sum(depth_values) / len(depth_values), 2)
        max_depth = max(depth_values)
        special_reads_file_content.write(f">{scaffold}:{start}-{end};{avg_depth};{max_depth}\n")
        special_reads_file_content.write(f"{sequence}\n")
        print_info(f"Outputted sequence for scaffold {scaffold} from {start} to {end} (Avg Depth: {avg_depth}, Max Depth: {max_depth})")

def execute_extract_special_sequences(bam_file_path:str, output_folder:str, depth_threshold: int = DEPTH_THRESHOLD, minimum_sequence_length: int = MINIMUM_SEQUENCE_LENGTH, threads: int = THREADS_DEFAULT):
    
    output_filename = os.path.join(output_folder, os.path.basename(bam_file_path).replace(FILE_ENDING_BAM, f"_special_reads_depth_gt_{DEPTH_THRESHOLD}.fasta"))
    
    print_info(f"Processing BAM file: {bam_file_path} with depth threshold: {depth_threshold}. Outputting to: {output_filename}")
    
    output_file_path = os.path.join(output_folder, output_filename)

    if os.path.exists(output_file_path):
        print_skipping(f"Output file {output_file_path} already exists!")
        return

    with pysam.AlignmentFile(bam_file_path, 'rb', threads=threads) as bamfile, open(output_file_path, 'w') as special_reads_file_content:

        current_sequence = ""
        current_start_position = None
        current_depth_values = []
        current_scaffold = None
        total_sequences = 0

        for pileup_column in bamfile.pileup():
            position = pileup_column.reference_pos + 1
            scaffold = bamfile.get_reference_name(pileup_column.reference_id)
            depth = pileup_column.nsegments

            # If we are at a new scaffold, or depth is below threshold, output the current sequence (if any)
            if scaffold != current_scaffold or depth < depth_threshold:
                if current_sequence:
                    write_fasta_entry(special_reads_file_content, current_scaffold, current_start_position, position - 1, current_sequence, current_depth_values, minimum_sequence_length)
                    total_sequences += 1                    
                    current_sequence = ""  # Reset sequence
                    current_start_position = None
                    current_depth_values = []

            if depth >= depth_threshold:
                # Start a new sequence if we're not currently building one
                if not current_sequence:
                    current_start_position = position
                    current_scaffold = scaffold

                # Count occurrences of each base at this position
                bases = [pileup_read.alignment.query_sequence[pileup_read.query_position]
                         for pileup_read in pileup_column.pileups if pileup_read.query_position is not None]
                if bases:
                    most_common_base, _ = Counter(bases).most_common(1)[0]
                    current_sequence += most_common_base
                    current_depth_values.append(depth)
            else:
                current_scaffold = scaffold  # Update the scaffold

        # If a sequence is still being built at the end, output it
        if current_sequence:
            write_fasta_entry(special_reads_file_content, current_scaffold, current_start_position, position - 1, current_sequence, current_depth_values, minimum_sequence_length)                
            total_sequences += 1

    print_success(f"Total sequences written for {bam_file_path}: {total_sequences}\n")

def write_unmapped_region(txtfile, scaffold: str, start: int, end: int, minimum_sequence_length: int=MINIMUM_SEQUENCE_LENGTH, maximum_sequence_length: int=MAXIMUM_SEQUENCE_LENGTH):
    """Writes an unmapped region entry to the output FASTA file."""
    if start is not None:
        
        if end - start + 1 < minimum_sequence_length or end - start + 1 > maximum_sequence_length:
            return

        txtfile.write(f"{scaffold}:{start}-{end}\n")
        print_info(f"Unmapped region found in {scaffold}: {start}-{end}")

def execute_extract_unmapped_regions(bam_file_path: str, output_folder: str, minimum_sequence_length: int = MINIMUM_SEQUENCE_LENGTH, maximum_sequence_length: int = MAXIMUM_SEQUENCE_LENGTH, threads: int = THREADS_DEFAULT):
    """Extracts regions of the reference genome with no coverage (depth = 0) from a BAM file."""
    output_filename = os.path.join(output_folder, os.path.basename(bam_file_path).replace(FILE_ENDING_BAM, "_unmapped_regions.txt"))
    
    if os.path.exists(output_filename):
        print_skipping(f"Output file {output_filename} already exists!")
        return

    print_info(f"Extracting unmapped regions from {bam_file_path}. Output: {output_filename}")

    with pysam.AlignmentFile(bam_file_path, 'rb', threads=threads) as bamfile, open(output_filename, 'w') as txtfile:
        total_regions = 0
        
        # Iterate over all reference scaffolds
        for scaffold, scaffold_length in zip(bamfile.references, bamfile.lengths):
            uncovered_start = None  # Start of an uncovered region

            # Convert covered positions into a set for fast lookup
            covered_positions = {pileup.reference_pos + 1 for pileup in bamfile.pileup(scaffold)}

            for position in range(1, scaffold_length + 1):  # 1-based positions
                if position not in covered_positions:
                    if uncovered_start is None:
                        uncovered_start = position  # Start of a new uncovered region
                else:
                    if uncovered_start is not None:
                        # End of an uncovered region found
                        write_unmapped_region(txtfile, scaffold, uncovered_start, position - 1, minimum_sequence_length, maximum_sequence_length)
                        total_regions += 1
                        uncovered_start = None

            # If the scaffold ends with an uncovered region, write it out
            if uncovered_start is not None:
                write_unmapped_region(txtfile, scaffold, uncovered_start, scaffold_length, minimum_sequence_length, maximum_sequence_length)
                total_regions += 1

    print_info(f"Total unmapped regions written for {bam_file_path}: {total_regions}\n")


def extract_special_sequences_for_species(species: str, depth_threshold: int = DEPTH_THRESHOLD):
    print_species_execution(f"Extracting special sequences for species {species} with depth threshold: {depth_threshold}")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.EXTRACT_SPECIAL_SEQUENCES,
        species,
        default=False
    ):
        print_skipping("Extracting special sequences is not enabled in the config. Skipping this step.")
        return

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_tuple in ref_genome_list:

        # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
        ref_genome_id = ref_genome_tuple[0]
        #ref_genome_path = ref_genome_tuple[1]

        print_info(f"Extracting special sequences for species {species} and reference genome {ref_genome_id}...")

        mapped_folder = get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)
        target_folder = get_folder_path_species_results_refgenome_special_sequences(species, ref_genome_id)

        bam_files = get_files_in_folder_matching_pattern(mapped_folder, f"*{FILE_ENDING_SORTED_BAM}")

        if len(bam_files) == 0:
            print_warning(f"No BAM files found in {mapped_folder} for species {species}.")
            return

        print_debug(f"Found {len(bam_files)} BAM files.")
        print_debug(f"BAM files: {bam_files}")

        for bam_file in bam_files:
            execute_extract_special_sequences(bam_file, target_folder, depth_threshold)
            execute_extract_unmapped_regions(bam_file, target_folder, 1000, 10000)
    
    print_info(f"Finished extracting special sequences for species {species}")

def all_species_extract_special_sequences(depth_threshold: int = DEPTH_THRESHOLD):
    print_step_execution("Extract special sequences for all species")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.PREPARE_REFERENCE_GENOME,
        default=False
    ):
        print_skipping("Extracting special sequences is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 
        extract_special_sequences_for_species(species, depth_threshold)

    print_success("Finished extracting special sequences for all species")
