
import os
import glob

# -----------------------------------------------------------------------------------------------
# Get all files in a folder matching a specific pattern (e.g., *.fastq.gz)
def get_files_in_folder_matching_pattern(folder: str, pattern: str) -> list:
    # Check if the folder exists
    if not os.path.exists(folder):
        raise Exception(f"Invalid folder: {folder}")
    # Read all files matching the pattern into a list
    files = glob.glob(os.path.join(folder, pattern))
    return files

# -----------------------------------------------------------------------------------------------
# Get all raw read files for a given species    
def get_read_files_for_species(species: str) -> list[str]:

    try:
        read_folder = os.path.join(species, "raw", "reads")
        read_files = get_files_in_folder_matching_pattern(read_folder, "*.fastq.gz")
    except Exception as e:  
        # Try looking in species folder directly as fallback.
        read_files = get_files_in_folder_matching_pattern(species, "*.fastq.gz")
        
    if len(read_files) == 0:
        raise Exception(f"No read files found for species {species}.")
        
    return read_files

# -----------------------------------------------------------------------------------------------
# Get only R1 raw read files for a given species
def get_r1_read_files_for_species(species: str) -> list[str]:
    files = get_read_files_for_species(species)
    r1_files = [f for f in files if "_R1" in os.path.basename(f)]
    if len(r1_files) == 0:
        raise Exception(f"No R1 read files found for species {species}.")
    return r1_files

# -----------------------------------------------------------------------------------------------
# Get sample IDs for a species based on raw read filenames
def get_sample_ids_for_species(species):
    files = get_r1_read_files_for_species(species)

    samples = []
    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        samples.append(filename)
    
    return samples

# -----------------------------------------------------------------------------------------------
# Get individual sample IDs for a given species based on raw read filenames
def get_individuals_for_species(species):

    samples = get_sample_ids_for_species(species)
    if len(samples) == 0:
        raise Exception(f"No raw reads found for species {species}.")
    individuals = set()
    for s in samples:
        individuals.add(get_individual_from_sample(s))
    return sorted(list(individuals))

# -----------------------------------------------------------------------------------------------
# Get reference genome files for a given species (supports .fna, .fasta, .fa)
def get_reference_genome_file_list_for_species(species: str) -> list[tuple[str, str]]:
    # Construct reference genome folder path
    species_folder = species
    ref_genome_folder = os.path.join(species, "raw", "ref_genome")
    try:
    # Collect all supported reference genome files
        reference_genome_files = get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fna")
        reference_genome_files += get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fasta")
        reference_genome_files += get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fa")
    except Exception as e:
        # Try looking in species folder directly as fallback.
        reference_genome_files = get_files_in_folder_matching_pattern(species_folder, f"*.fna")
        reference_genome_files += get_files_in_folder_matching_pattern(species_folder, f"*.fasta")
        reference_genome_files += get_files_in_folder_matching_pattern(species_folder, f"*.fa")

    if len(reference_genome_files) == 0:
        raise Exception(f"No reference genome found for species {species}.")
        
    # Return as list of tuples: (filename without extension, full path)
    reference_genome_files_with_filename = [(os.path.splitext(os.path.basename(f))[0], f) for f in reference_genome_files]
    return reference_genome_files_with_filename

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
# Get only reference genome file paths for a species
def get_ref_genomes_for_species(species):
    refs = get_reference_genome_file_list_for_species(species)
    return [ref[1] for ref in refs]
