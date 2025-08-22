import os
import subprocess
import glob
from typing import Optional

from common.common_constants import *
from common.common_logging import *
from common.common_config import *
from common.common_folder_functions import *
from common.common_config_enumerations import *

#####################
# Helpers
#####################

def is_species(species: str) -> bool:
    return species in config.get_species_list()

def get_species_name(species_id):

    if not is_species(species_id):
        raise ValueError(f"Invalid species ID: {species_id}")
   
    return config.get_species_setting(species_id).name


def is_sam_file_sorted(sam_file: str) -> bool:
    """
    Checks if a SAM file is sorted by coordinate.

    :param sam_file: The path to the SAM file
    :return: True if the SAM file is sorted by coordinate, False otherwise
    """
    try:
        # Check the header for sorting status
        result = subprocess.run(
            [PROGRAM_PATH_SAMTOOLS, PROGRAM_PATH_SAMTOOLS_VIEW, '-H', sam_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        for line in result.stdout.splitlines():
            if line.startswith('@HD') and 'SO:coordinate' in line:
                # The SAM file is sorted by coordinate
                return True
        # The SAM file is not sorted by coordinate
        return False
    except subprocess.CalledProcessError as e:
        print_error(f"Error reading header: {e}")
        # Return False in case of an error
        return False

def index_fasta_file_samtools(fasta_file_path: str):
   # Index the newly created gzipped consensus sequence using Samtools faidx.

    print_info(f"SAMTOOLS indexing reference genome {fasta_file_path} ...")

   # check if the index file already exists
    index_file = f"{fasta_file_path}.fai"
    if os.path.exists(index_file):
        print_skipping(f"Index file {index_file} already exists.")
        return
   
    command_samtools = [
        PROGRAM_PATH_SAMTOOLS, 
        PROGRAM_PATH_SAMTOOLS_FAIDX, 
        fasta_file_path
        ]

    run_command(command_samtools, 
                description=f"Indexing fasta file {fasta_file_path}"
    )

def index_fasta_file_bwa(fasta_file_path: str):
    print_info(f"BWA Indexing reference genome {fasta_file_path} ...")

    # Check if index file already exists
    index_file = f"{fasta_file_path}.bwt"
    if os.path.exists(index_file):
        print_skipping(f"Index file {index_file} already exists.")
        return 

    command_bwa = [PROGRAM_PATH_BWA, PROGRAM_PATH_BWA_INDEX, fasta_file_path ]

    run_command(command_bwa, 
            description=f"Indexing fasta file {fasta_file_path}"
    )

def is_fasta_file(file_name: str) -> bool:
    return file_name.endswith(FILE_ENDING_FNA) or file_name.endswith(FILE_ENDING_FA) or file_name.endswith(FILE_ENDING_FASTA)

def is_fasta_gz_file(file_name: str) -> bool:
    return file_name.endswith(FILE_ENDING_FNA_GZ) or file_name.endswith(FILE_ENDING_FA_GZ) or file_name.endswith(FILE_ENDING_FASTQ_GZ) 

def call_r_script(script_path: str, *args):
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"R script not found: {script_path}")

    command = ["Rscript", script_path] + list(args)

    run_command(command, description=f"R script: {script_path}", throw_error=False)

def get_adapter_sequence(species: str) -> tuple[str,str]:
    """
    Get the adapter sequence for a given species. If not found, use the global adapter sequence.

    :param species: The species name
    :return: The adapter sequence for R1 and R2
    """
    # Check if the species is valid
    if not is_species(species):
        raise ValueError(f"Invalid species: {species}")
    
    adapter_removal_config = config.get_pipeline_settings(PipelineStages.RAW_READS_PROCESSING, RawReadsProcessingSteps.ADAPTER_REMOVAL, species)    
   
    print(adapter_removal_config)

    adapter_sequence_r1 = adapter_removal_config['settings']['adapters_sequences']["r1"]
    adapter_sequence_r2 = adapter_removal_config['settings']['adapters_sequences']["r2"]

    if adapter_sequence_r1 is None and adapter_sequence_r2 is None:
        raise ValueError(f"Adapter sequence not found for species: {species}")

    return adapter_sequence_r1, adapter_sequence_r2

def get_filename_from_path(file_path: str) -> str:
    return os.path.basename(file_path)

def get_filename_from_path_without_extension(file_path: str) -> str:
    return os.path.splitext(get_filename_from_path(file_path))[0]

def get_r_script(r_script_name: str, processing_folder: str) -> str:
    r_script_path = os.path.join(get_folder_path_scripts_plots(processing_folder), r_script_name)

    if not os.path.exists(r_script_path):
        raise Exception(f"Invalid R script: {r_script_path}")   
    
    return r_script_path

def get_raw_reads_list_of_species(species: str) -> list:
     
    if not is_species_folder(species):
        raise Exception(f"Invalid species folder: {species}")
    
    raw_reads_folder = get_folder_path_species_raw_reads(species)

    return get_files_in_folder_matching_pattern(raw_reads_folder, f"*{FILE_ENDING_FASTQ_GZ}")
   

def get_files_in_folder_matching_pattern(folder: str, pattern: str) -> list:
     
    if not os.path.exists(folder):
        raise Exception(f"Invalid folder: {folder}")
    
    #read all reads from folder into list
    files = glob.glob(os.path.join(folder, pattern))

    return files

def get_raw_paired_reads_list_of_species(species: str) -> list:

    if not is_species_folder(species):
        raise Exception(f"Invalid species folder: {species}")

    folder_path = get_folder_path_species_raw_reads(species)

    find_command = f'find {folder_path} -type f \\( -name "{FILE_PATTERN_R1_FASTQ_GZ}" -o -name "{FILE_PATTERN_R2_FASTQ_GZ}" \\) ! -path "{FILE_PATTERN_UNDETERMINED_FOLDER}"'
    
    # Sort the list of files by sample ID
    sort_command = 'sort'
    
    # Use awk to format the list of files as a comma-separated list
    awk_command = 'awk \'NR%2{printf "%s,", $0} NR%2==0{print $0}\''
    
    # Combine the commands into a single command string
    full_command = f'{find_command} | {sort_command} | {awk_command}'

    print_debug(f"Running command: {full_command}")
    
    try:
        #print_info(f"Running command: {full_command}")
        command_result = subprocess.run(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

    except subprocess.CalledProcessError as e:
        raise Exception(f"An error occurred while running the command: {e}")
        
    #convert output into list of files
    lines = command_result.stdout.splitlines()

    # Split each line by the comma, strip the paths to remove any extra spaces or newlines, and make them absolute paths
    file_paths = [
        [os.path.abspath(os.path.join(folder_path, paths[0].strip())),
        os.path.abspath(os.path.join(folder_path, paths[1].strip()))]
        for line in lines
        for paths in [line.strip().split(',')]  # Split the line by the comma
    ]

    return file_paths

def run_command(command: list, description: str = "", cwd: Optional[str] = None, throw_error: bool = True) -> str:
    """
    Run a shell command and return its stdout output.
    Logs stderr normally unless result.returncode != 0 and throw_error=True.
    """
    print_debug("Entering run_command function")

    command_program = command[0] if command else "Unknown"

    print_info(f"Running: {description or command_program}")
    print_debug(f"Command program: {command_program}")
    print_debug(f"Command: {' '.join(command)}")

    if cwd:
        print_debug(f"Working directory: {cwd}")

    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=cwd
    )

    # Print stdout
    if result.stdout:
        for line in result.stdout.splitlines():
            print_info(f"[{command_program}] {line.strip()}")

    # Print stderr
    if result.stderr:
        for line in result.stderr.splitlines():
            if result.returncode == 0 or not throw_error:
                print_info(f"[{command_program}] {line.strip()}")  # Treat as log
            else:
                print_error(f"[{command_program}] {line.strip()}")  # Real error

    # Raise exception only if returncode != 0 and throw_error is True
    if result.returncode != 0:
        print_error(f"{command_program} exited with return code {result.returncode}")
        if throw_error:
            raise subprocess.CalledProcessError(result.returncode, command, output=result.stdout, stderr=result.stderr)

    print_info(f"{command_program} completed with return code {result.returncode}")
    return result.stdout.strip()

