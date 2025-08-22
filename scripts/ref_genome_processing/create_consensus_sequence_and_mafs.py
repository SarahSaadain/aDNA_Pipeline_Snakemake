import os
import subprocess

from multiprocessing import Pool, cpu_count
from common_aDNA_scripts import *
import ref_genome_processing.common_ref_genome_processing_helpers as common_rgp

from ref_genome_processing.map_aDNA_to_refgenome import execute_bwa_map_aDNA_to_refgenome

def execute_angsd_create_consensus_sequence(sorted_bam_file: str, reference_genome_path: str, output_dir: str):

    pid = os.getpid() # Get current process ID for logging
    print_info(f"[PID {pid}] Executing angsd to create snps and consensus sequence for {sorted_bam_file}")

    if not os.path.exists(sorted_bam_file):
        print_warning(f"[PID {pid}] Input file {sorted_bam_file} does not exist!")
    
    if not os.path.exists(output_dir):
        print_error(f"[PID {pid}] Output directory {output_dir} does not exist!")

    base_name = get_filename_from_path_without_extension(sorted_bam_file)

    # Construct the output file path for the consensus FASTA sequence.
    consensus_out_file_path = os.path.join(output_dir, f"{base_name}_consensus")
    snp_out_file_path = os.path.join(output_dir, f"{base_name}_snp")

    print_debug(f"[PID {pid}] Output file path consensus: {consensus_out_file_path}")
    print_debug(f"[PID {pid}] Output file path snp: {snp_out_file_path}")
    
    print_info(f"[PID {pid}] Creating snps and consensus sequence of {sorted_bam_file}...")

    command_angsd_consensus = [
        PROGRAM_PATH_ANGSD,
        "-out", consensus_out_file_path,
        "-i", sorted_bam_file,
        "-ref", reference_genome_path,
        "-doFasta", "2",            # Use most common base
        "-doCounts", "1",           # Needed for doFasta
        "-minMapQ", "10",
        "-minQ", "0",               # Keep rescaled bases
        "-remove_bads", "1",
        "-baq", "1",
        "-C", "50"
    ]

    command_angsd_snp = [
        PROGRAM_PATH_ANGSD,
        "-out", snp_out_file_path,
        "-i", sorted_bam_file,
        "-ref", reference_genome_path,
        "-doMaf", "1",
        "-doMajorMinor", "1",
        "-GL", "2",
        "-SNP_pval", "1e-6",
        #"-dosnpstat", "1",
        "-minMapQ", "30",
        "-minQ", "20",                 # Keep strict filtering for SNPs
        "-remove_bads", "1",
        "-baq", "1",
        "-C", "50"
    ]


    try:

        if os.path.exists(snp_out_file_path + ".mafs.gz"):
            print_skipping(f"[PID {pid}] SNPs for {base_name} already exists.")
        else:

            run_command(
                command_angsd_snp, 
                description=f"[PID {pid}] Creating SNPs for {sorted_bam_file}"
            )

        # Check if the gzipped consensus sequence already exists to avoid re-processing.
        if os.path.exists(consensus_out_file_path + FILE_ENDING_FA_GZ):
            print_skipping(f"[PID {pid}] Consensus sequence for {base_name} already exists.")
        else:

            run_command(
                command_angsd_consensus, 
                description=f"[PID {pid}] Creating consensus sequence for {sorted_bam_file}"
            )
    
    except subprocess.CalledProcessError as e:
        print_error(f"[PID {pid}] Failed to create snps and consensus sequence for {sorted_bam_file} : {e}")
        return


def create_consensus_sequence_for_species_and_ref(species: str, ref_genome_tuple: tuple[str, str]):

    ref_genome_id = ref_genome_tuple[0]
    ref_genome_path = ref_genome_tuple[1]

    print_info(f"Creating consensus sequence for species {species} and ref genome {ref_genome_id}...")

    # Get the folder containing mapped aDNA reads for the given species and reference genome.
    aDNA_reads_folder = get_folder_path_species_processed_refgenome_mapped(species, ref_genome_id)

    # Get a list of all sorted BAM files (mapped reads) in the specified folder.
    sorted_bam_paths = get_files_in_folder_matching_pattern(aDNA_reads_folder, f"*{FILE_ENDING_SORTED_BAM}")

    # If no mapped reads are found, log a warning and exit the function.
    if len(sorted_bam_paths) == 0:
        print_warning(f"No mapped reads found for species {species} with ref genome {ref_genome_id}. Skipping consensus sequence creation.")
        return
    
    print_debug(f"Found {len(sorted_bam_paths)} mapped reads for species {species} with ref genome {ref_genome_id}.")
    print_debug(f"Mapped reads: {sorted_bam_paths}")

    # Determine the output folder for consensus sequences.
    output_folder_angsd = get_folder_path_species_processed_refgenome_consensus_sequences(species, ref_genome_id)

    files_in_folder = get_files_in_folder_matching_pattern(output_folder_angsd, f"*")

    # Determine the number of processes to use for parallel execution.
    # It takes the minimum of a predefined default (THREADS_DEFAULT) and the
    # actual number of CPU cores available on the system, ensuring efficient
    # resource utilization without overloading the system.
    num_processes = min(THREADS_DEFAULT, cpu_count(), len(sorted_bam_paths))
    print_debug(f"Using {num_processes} processes for parallel execution.")

    args_for_pool = []

    # Initialize a multiprocessing Pool. This creates a pool of worker processes
    # that can execute tasks concurrently. The 'processes' argument specifies
    # the maximum number of worker processes to use.
    with Pool(processes=num_processes) as pool:
        # Prepare the arguments for each call to the _process_single_mapped_aDNA_file helper function.
        # Each tuple (mapped_aDNA_read_file_path, species, ref_genome_id) represents
        # one set of arguments for a single task.
        for sorted_bam_file in sorted_bam_paths:
            
            sorted_rescaled_bam_file_path = common_rgp.get_rescaled_bam_path_for_sorted_bam_path(sorted_bam_file)

            # use the rescaled file if it exists
            if os.path.exists(sorted_rescaled_bam_file_path):
                print_info(f"Using rescaled file {sorted_rescaled_bam_file_path} for consensus sequence creation.")
                sorted_bam_file = sorted_rescaled_bam_file_path

            args_for_pool.append((sorted_bam_file, ref_genome_path, output_folder_angsd))

        #args_for_pool = [(mapped_aDNA_read_file_path, ref_genome_path, output_folder_angsd) for mapped_aDNA_read_file_path in list_of_mapped_aDNA_files]

        # Use pool.starmap to apply the _process_single_mapped_aDNA_file function to each
        # set of arguments in args_for_pool. starmap is suitable when the target
        # function expects multiple arguments, which are provided as a tuple.
        # This call blocks until all tasks in the pool have completed.
        pool.starmap(execute_angsd_create_consensus_sequence, args_for_pool)

    print_info(f"Consensus sequence creation for species {species} and ref genome {ref_genome_id} complete.")

def create_consensus_sequence_for_species(species):

    print_species_execution(f"Creating consensus sequence for species {species}")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.CREATE_CONSENSUS_SEQUENCE,
        species,
        default=False
    ):
        print_skipping(f"Creating consensus sequence for species {species} is disabled in the config. Skipping this species.")
        return

    try:
        ref_genome_list = common_rgp.get_reference_genome_file_list_for_species(species)
    except Exception as e:
        print_error(f"Failed to get reference genome files for species {species}: {e}")
        return

    for ref_genome_tuple in ref_genome_list:

        # add x/x message
        print_info(f"[{ref_genome_list.index(ref_genome_tuple) + 1}/{len(ref_genome_list)}] Creating consensus sequence for species {species} and ref genome {ref_genome_tuple[0]}...")
    
        create_consensus_sequence_for_species_and_ref(species, ref_genome_tuple)      

    print_info(f"Consensus sequence of {species} created and mapped successfully.")

def all_species_create_consensus_sequence_and_mafs():

    print_step_execution("Creating and mapping consensus sequence for all species")

    if not config.is_process_step_enabled(
        PipelineStages.REFERENCE_GENOME_PROCESSING, 
        ReferenceGenomeProcessingSteps.CREATE_CONSENSUS_SEQUENCE
    ):
        print_skipping("Creating consensus sequences is not enabled in the config. Skipping this step.")
        return

    for species in FOLDER_SPECIES: 

        create_consensus_sequence_for_species(species)

    print_info("Creating and mapping consensus sequence for all species complete")