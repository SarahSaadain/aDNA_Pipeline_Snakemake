import os
import subprocess

from multiprocessing import Pool, cpu_count
from common_aDNA_scripts import *
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

from ref_genome_processing.map_aDNA_to_refgenome import execute_bwa_map_aDNA_to_refgenome

def execute_angsd_create_and_map_consensus_sequence(sorted_bam_file: str, output_dir: str):

    pid = os.getpid() # Get current process ID for logging
    print_info(f"[PID {pid}] Executing angsd to create consensus sequence for {sorted_bam_file}")

    if not os.path.exists(sorted_bam_file):
        print_warning(f"[PID {pid}] Input file {sorted_bam_file} does not exist!")
    
    if not os.path.exists(output_dir):
        print_error(f"[PID {pid}] Output directory {output_dir} does not exist!")

    base_name = get_filename_from_path_without_extension(sorted_bam_file)

    # Construct the output file path for the consensus FASTA sequence.
    out_file_path = os.path.join(output_dir, f"{base_name}_consensus{FILE_ENDING_FASTA}")

    print_debug(f"[PID {pid}] Output file path: {out_file_path}")
    
    # Check if the gzipped consensus sequence already exists to avoid re-processing.
    if os.path.exists(out_file_path + FILE_ENDING_FA_GZ):
        print_info(f"[PID {pid}] Consensus sequence for {base_name} already exists. Skipping.")
        return

    print_info(f"[PID {pid}] Creating consensus sequence of {sorted_bam_file}...")

    command_angsd = [
        PROGRAM_PATH_ANGSD, 
        "-out", out_file_path, 
        "-i", sorted_bam_file, 
        "-doFasta", "2", 
        "-doCounts", "1"
        "-doMajorMinor", "1"
        "-doMaf", "1",
        "-GL", "1",
        "-SNP_pval", "1e-6"
    ]

    try:

        run_command(command_angsd, 
                    description=f"[PID {pid}] Creating consensus sequence for {sorted_bam_file}",
                    throw_error=True
        )
    
    except subprocess.CalledProcessError as e:
        print_error(f"[PID {pid}] Failed to create consensus sequence for {sorted_bam_file} : {e}")
        return

    # Index the newly created gzipped consensus sequence using Samtools faidx.
    print_info(f"[PID {pid}] Indexing consensus sequence {out_file_path + FILE_ENDING_FA_GZ}...")

    command_samtools = [
        PROGRAM_PATH_SAMTOOLS, PROGRAM_PATH_SAMTOOLS_FAIDX, 
        "-i", out_file_path+FILE_ENDING_FA_GZ
        ]
    
    try:
        run_command(command_samtools, 
                    description=f"[PID {pid}] Indexing consensus sequence for {sorted_bam_file}",  
                    throw_error=True
        )

        print_info(f"[PID {pid}] Consensus sequence {out_file_path + FILE_ENDING_FA_GZ} indexed successfully.")
    
    except subprocess.CalledProcessError as e:
        print_error(f"[PID {pid}] Failed to index consensus sequence {out_file_path}: {e}")


def create_consensus_sequence_for_species(species: str, ref_genome_id: str):

    print_info(f"Creating consensus sequence for species {species} and ref genome {ref_genome_id}...")

    # Get the folder containing mapped aDNA reads for the given species and reference genome.
    aDNA_reads_folder = get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)

    # Get a list of all sorted BAM files (mapped reads) in the specified folder.
    list_of_mapped_aDNA_files = get_files_in_folder_matching_pattern(aDNA_reads_folder, f"*{FILE_ENDING_SORTED_BAM}")

    # If no mapped reads are found, log a warning and exit the function.
    if len(list_of_mapped_aDNA_files) == 0:
        print_warning(f"No mapped reads found for species {species} with ref genome {ref_genome_id}. Skipping consensus sequence creation.")
        return
    
    print_debug(f"Found {len(list_of_mapped_aDNA_files)} mapped reads for species {species} with ref genome {ref_genome_id}.")
    print_debug(f"Mapped reads: {list_of_mapped_aDNA_files}")

    # Determine the output folder for consensus sequences.
    output_folder_consensus_seq = get_folder_path_species_processed_refgenome_mtdna_consensus_sequences(species, ref_genome_id)

    # Determine the number of processes to use for parallel execution.
    # It takes the minimum of a predefined default (THREADS_DEFAULT) and the
    # actual number of CPU cores available on the system, ensuring efficient
    # resource utilization without overloading the system.
    num_processes = min(THREADS_DEFAULT, cpu_count())
    print_debug(f"Using {num_processes} processes for parallel execution.")

    # Initialize a multiprocessing Pool. This creates a pool of worker processes
    # that can execute tasks concurrently. The 'processes' argument specifies
    # the maximum number of worker processes to use.
    with Pool(processes=num_processes) as pool:
        # Prepare the arguments for each call to the _process_single_mapped_aDNA_file helper function.
        # Each tuple (mapped_aDNA_read_file_path, species, ref_genome_id) represents
        # one set of arguments for a single task.
        args_for_pool = [(mapped_aDNA_read_file_path, output_folder_consensus_seq) for mapped_aDNA_read_file_path in list_of_mapped_aDNA_files]

        # Use pool.starmap to apply the _process_single_mapped_aDNA_file function to each
        # set of arguments in args_for_pool. starmap is suitable when the target
        # function expects multiple arguments, which are provided as a tuple.
        # This call blocks until all tasks in the pool have completed.
        pool.starmap(execute_angsd_create_and_map_consensus_sequence, args_for_pool)

    print_info(f"Consensus sequence creation for species {species} and ref genome {ref_genome_id} complete.")

def map_consensus_sequence_for_species(species: str, ref_genome_tuple: tuple):
    print_info(f"Mapping aDNA consensus sequence to reference genome for species {species} ...")

     # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
    ref_genome_id = ref_genome_tuple[0]
    ref_genome_path = ref_genome_tuple[1]

    #get reads
    read_folder = get_folder_path_species_processed_refgenome_consensus_sequences(species, ref_genome_id)
    list_of_read_files = get_files_in_folder_matching_pattern(read_folder, f"*{FILE_ENDING_FASTQ_GZ}")

    if len(list_of_read_files) == 0:
        print_warning(f"No reads found for species {species}. Skipping.")
        return
    
    print_debug(f"Found {len(list_of_read_files)} reads for species {species}.")
    print_debug(f"Reads: {list_of_read_files}")
    
    output_folder = get_folder_path_species_processed_refgenome_mtdna_consensus_sequences_mapped(species, ref_genome_id)

    for read_file_path in list_of_read_files:

        print_info(f"Mapping {read_file_path} to reference genome {ref_genome_path} ...")

        read_name = os.path.splitext(os.path.basename(read_file_path))[0]
        output_file_path = os.path.join(output_folder, f"{read_name}_{ref_genome_id}{FILE_ENDING_SAM}")

        execute_bwa_map_aDNA_to_refgenome(read_file_path, ref_genome_path, output_file_path, THREADS_DEFAULT)

    print_info(f"Mapping aDNA consensus sequence to reference genome for species {species} complete")

def create_and_map_consensus_sequence_for_species(species):

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_tuple in ref_genome_list:

        # ref_genome is a tuple of (ref_genome_name without extension, ref_genome_path)
        ref_genome_id = ref_genome_tuple[0]
        #ref_genome_path = ref_genome_tuple[1]
    
        #create_consensus_sequence_for_species(species, ref_genome_id)
        map_consensus_sequence_for_species(species, ref_genome_tuple)        

    print_info(f"Consensus sequence of {species} created and mapped successfully.")

def all_species_create_and_map_consensus_sequence():
    print_step_execution("Creating and mapping consensus sequence for all species")

    for species in FOLDER_SPECIES: 
        create_and_map_consensus_sequence_for_species(species)

    print_info("Creating and mapping consensus sequence for all species complete")