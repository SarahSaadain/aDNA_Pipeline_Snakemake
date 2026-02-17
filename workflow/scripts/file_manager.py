
import os
import glob
import re

# -----------------------------------------------------------------------------------------------
# Get all files in a folder matching a specific pattern (e.g., *.fastq.gz)
def get_files_in_folder_matching_pattern(folder: str, pattern: str) -> list:
    # Check if the folder exists
    if not os.path.exists(folder):
        logger.error(f"Invalid folder: {folder}")
        raise Exception(f"Invalid folder: {folder}")
    # Read all files matching the pattern into a list
    files = glob.glob(os.path.join(folder, pattern))
    return files

# -----------------------------------------------------------------------------------------------
# Get all raw read files for a given species    
def get_read_files_for_species(species: str) -> list[str]:

    read_folder = os.path.join(species, "raw", "reads")

    try:   
        logger.debug(f"Looking for read files in {read_folder} for species {species}.")
        read_files = get_files_in_folder_matching_pattern(read_folder, "*.fastq.gz")
    except Exception as e:  
        # Try looking in species folder directly as fallback.
        logger.debug(f"Read folder not found for species {species}. Trying species folder directly.")
        read_files = get_files_in_folder_matching_pattern(species, "*.fastq.gz")
        
    if len(read_files) == 0:
        logger.error(f"No read files found for species {species}.")
        raise Exception(f"No read files found for species {species}.")
    
    logger.debug(f"Read files for species {species}: {read_files}")
        
    return read_files

# -----------------------------------------------------------------------------------------------
# Get only R1 raw read files for a given species
def get_r1_read_files_for_species(species: str) -> list[str]:
    files = get_read_files_for_species(species)
    r1_files = [f for f in files if "_R1" in os.path.basename(f)]
    if len(r1_files) == 0:
        logger.error(f"No R1 read files found for species {species}.")
        logger.error(f"Available read files for species {species}: {files}")
        raise Exception(f"No R1 read files found for species {species}.")
    
    logger.debug(f"R1 read files for species {species}: {r1_files}")
    
    return r1_files

# -----------------------------------------------------------------------------------------------
# Get sample IDs for a species based on raw read filenames
def get_sample_ids_for_species(species):
    files = get_r1_read_files_for_species(species)

    samples = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        samples.append(filename)
    
    logger.debug(f"Sample IDs for species {species}: {samples}")
    
    return samples

def get_raw_reads_for_sample(species, sample):

    read_files = get_read_files_for_species(species) 
    
    # turn read paths into file names only
    read_files = [os.path.basename(f) for f in read_files]

    reads_dir = f"{species}/raw/reads"
 
    # R1
    base_r1 = f"{sample}_R1"
    candidates_r1 = [f for f in read_files
                     if re.match(base_r1 + r"(\S*)?\.fastq\.gz", f)]
    
    if not candidates_r1:
        logger.error(f"No R1 found for {sample}. Expected pattern: {base_r1}*.fastq.gz in {reads_dir}. Found files: {read_files}")
        raise FileNotFoundError(f"No R1 found for {sample}. Expected pattern: {base_r1}*.fastq.gz in {reads_dir}. Found files: {read_files}")

    r1 = os.path.join(reads_dir, sorted(candidates_r1)[0])
 
    # R2
    base_r2 = f"{sample}_R2"
    candidates_r2 = [f for f in read_files
                     if re.match(base_r2 + r"(\S*)?\.fastq\.gz", f)]
 
    if candidates_r2:
        r2 = os.path.join(reads_dir, sorted(candidates_r2)[0])
        return [r1, r2]  # Paired-end
    else:
        return [r1]      # Single-end

# -----------------------------------------------------------------------------------------------
# Get individual sample IDs for a given species based on raw read filenames
def get_individuals_for_species(species):

    samples = get_sample_ids_for_species(species)
    if len(samples) == 0:
        logger.error(f"No samples found for species {species}.")
        raise Exception(f"No samples found for species {species}.")
    individuals = set()
    for s in samples:
        individuals.add(get_individual_from_sample(s))

    logger.debug(f"Individuals for species {species}: {individuals}")

    return sorted(list(individuals))

# -----------------------------------------------------------------------------------------------
# Get reference files for a given species (supports .fna, .fasta, .fa)
def get_reference_file_list_for_species(species: str) -> list[tuple[str, str]]:
    # Construct reference folder path
    species_folder = species
    reference_folder = os.path.join(species, "raw", "ref")
    try:
        # Collect all supported reference files
        logger.debug(f"Looking for reference files in {reference_folder} for species {species}.")

        reference_files = get_files_in_folder_matching_pattern(reference_folder, "*.fna")
        reference_files += get_files_in_folder_matching_pattern(reference_folder, "*.fasta")
        reference_files += get_files_in_folder_matching_pattern(reference_folder, "*.fa")
    except Exception as e:
        # Try looking in species folder directly as fallback.
        logger.debug(f"Reference folder not found for species {species}. Trying species folder directly.")

        reference_files = get_files_in_folder_matching_pattern(species_folder, "*.fna")
        reference_files += get_files_in_folder_matching_pattern(species_folder, "*.fasta")
        reference_files += get_files_in_folder_matching_pattern(species_folder, "*.fa")

    if len(reference_files) == 0:
        raise Exception(f"No reference found for species {species}.")
        
    # Return as list of tuples: (filename without extension, full path)
    reference_files_with_filename = [(os.path.splitext(os.path.basename(f))[0].replace('.', '_'), f) for f in reference_files]

    logger.debug(f"Reference files for species {species}: {reference_files_with_filename}")

    return reference_files_with_filename

# -----------------------------------------------------------------------------------------------
# Extract individual ID from a given file path or sample name
def get_individual_from_filepath(filepath):
    basename = os.path.basename(filepath)
    return get_individual_from_sample(basename)

# -----------------------------------------------------------------------------------------------
# Extract individual ID from a sample name
def get_individual_from_sample(sample):
    return sample.split("_")[0]

# -----------------------------------------------------------------------------------------------
# Get only reference file paths for a species
def get_references_paths_for_species(species):
    refs = get_reference_file_list_for_species(species)
    return [ref[1] for ref in refs]

# -----------------------------------------------------------------------------------------------
# Get only reference IDs for a species
def get_references_ids_for_species(species):
    refs = get_reference_file_list_for_species(species)
    return [ref[0] for ref in refs]

# -----------------------------------------------------------------------------------------------
# Get sample IDs for a specific individual within a species
def get_samples_for_species_individual(species, individual):
    samples = get_sample_ids_for_species(species)
    samples_of_individual = [f for f in samples if f.startswith(f"{individual}")]

    logger.debug(f"Samples for individual {individual} in species {species}: {samples_of_individual}")
    
    return samples_of_individual