import os
import glob
import logging

def get_files_in_folder_matching_pattern(folder: str, pattern: str) -> list:
     
    if not os.path.exists(folder):
        raise Exception(f"Invalid folder: {folder}")
    
    #read all reads from folder into list
    files = glob.glob(os.path.join(folder, pattern))

    return files


def get_reference_genome_file_list_for_species(species: str) -> list[tuple[str, str]]:
 
    # get ref genome
    ref_genome_folder = os.path.join(species, "raw", "ref_genome")

    # add fna files to reference genome list
    reference_genome_files = get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fna")
    # add fasta files to reference genome list
    reference_genome_files += get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fasta")
    # add fa files to reference genome list
    reference_genome_files += get_files_in_folder_matching_pattern(ref_genome_folder, f"*.fa")

    if len(reference_genome_files) == 0:
        raise Exception(f"No reference genome found for species {species}.")
    
    #common.print_debug(f"Found {len(reference_genome_files)} reference genome files for species {species}.")
    #common.print_debug(f"Reference genome files: {reference_genome_files}")

    # return as tuple of (filename without extension, filepath)
    reference_genome_files_with_filename = [(os.path.splitext(os.path.basename(f))[0], f) for f in reference_genome_files]

    return reference_genome_files_with_filename

# Function to get individuals for a given species
def get_individuals_for_species(species):

    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")

    if len(files) == 0:
        raise Exception(f"No raw reads found for species {species}.")

    individuals = set()
    for f in files:
        basename = os.path.basename(f)
        # take everything before the first underscore as individual ID
        individual = basename.split("_")[0]
        individuals.add(individual)
    return sorted(list(individuals))

def get_ref_genomes_for_species(species):
    #only return file paths
    refs = get_reference_genome_file_list_for_species(species)
    return [ref[1] for ref in refs]

def get_sample_ids_for_species(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")

    samples = []

    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
        samples.append(filename)

    if samples:
        print(f"Found {len(samples)} samples for species {species}.")
        print (f"Samples: {samples}")
    
    return samples

def get_input_multiqc_raw(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")

    all_inputs = []

    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','')
            

        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "raw", f"{filename}_fastqc.html"))

        if os.path.exists(raw_file.replace("_R1_", "_R2_")):
            all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "raw", f"{filename.replace("_R1_", "_R2_")}_fastqc.html"))
            

    return all_inputs

def get_input_reads_processing_raw(species):
    

    all_inputs = []

    for sample in get_sample_ids_for_species(species):
        #{species}/processed/qualitycontrol/statistics/{sample}_raw_reads.count
        all_inputs.append(os.path.join(species, "processed", "qualitycontrol", "statistics", f"{sample}_all_steps_count.csv"))
       
    return all_inputs

def get_input_reads_processing_trimmed(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")

    all_inputs = []

    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
            
        #{species}/processed/qualitycontrol/statistics/{sample}_trimmed_reads.count
        all_inputs.append(os.path.join(species, "processed", "qualitycontrol", "statistics", f"{filename}_trimmed_reads.count"))
       
    return all_inputs

def get_input_reads_processing_quality_filtered(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")

    all_inputs = []

    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]
            
        #{species}/processed/qualitycontrol/statistics/{sample}_quality_filtered_reads.count
        all_inputs.append(os.path.join(species, "processed", "qualitycontrol", "statistics", f"{filename}_quality_filtered_reads.count"))
       
    return all_inputs


def get_input_multiqc_trimmed(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")

    all_inputs = []

    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]

        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "quality_filtered", f"{filename}_quality_filtered_fastqc.html"))

    return all_inputs

def get_input_multiqc_quality_filtered(species):
    files = get_files_in_folder_matching_pattern(os.path.join(species, "raw", "reads"), "*R1*.fastq.gz")

    all_inputs = []

    for raw_file in files:
        filename = os.path.basename(raw_file).replace('.fastq.gz','').split("_R1")[0]

        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "quality_filtered", f"{filename}_quality_filtered_fastqc.html"))

    return all_inputs

def get_input_multiqc_merged(species):
    individuals = get_individuals_for_species(species)

    all_inputs = []

    for individual in individuals:
        
        all_inputs.append(os.path.join(species, "results", "qualitycontrol", "fastqc", "merged", f"{individual}_fastqc.html"))

    return all_inputs


def get_all_inputs(wildcards):

    # Generate inputs for rule all
    all_inputs = []

    for species in config["species"]:
        species_folder = species

        #all_inputs += get_input_multiqc_raw(species)
        #all_inputs += get_input_multiqc_trimmed(species)
        #all_inputs += get_input_multiqc_quality_filtered(species)
        #all_inputs += get_input_multiqc_merged(species)

        # Add QC report
        all_inputs.append(os.path.join(species_folder, "results","qualitycontrol", f"quality_check_report_{species}.html"))


        #all_inputs += get_input_reads_processing_raw(species)
        #all_inputs += get_input_reads_processing_trimmed(species)
        #all_inputs += get_input_reads_processing_quality_filtered(species)

        #{species}/processed/qualitycontrol/statistics/{species}_all_steps_count.csv
        #all_inputs.append(os.path.join(species_folder, "results", "qualitycontrol", "statistics", f"{species}_reads_processing.csv"))

        for sample in get_sample_ids_for_species(species):
            #"{species}/processed/qualitycontrol/ecmsd/{sample}/mapping/Mito_summary.txt"
            all_inputs.append(os.path.join(species_folder, "results", "qualitycontrol", "ecmsd", sample, "mapping", "Mito_summary.txt"))

        all_inputs.append(os.path.join(species_folder, "results", "plots", f"{species}_read_count_comparison.png"))
        all_inputs.append(os.path.join(species_folder, "results", "plots", f"{species}_read_count_comparison_by_individual.png"))

        individuals = get_individuals_for_species(species)
        for ind in individuals:
            #all_inputs.append(os.path.join(species_folder, "processed", "merged", f"{ind}.fastq.gz"))

            #{species}/results/qualitycontrol/fastqc/merged/{individual}_merged_fastqc.html"
            try:
                ref_genome_list = get_reference_genome_file_list_for_species(species)
                
                for ref_genome_tuple in ref_genome_list:

                    ref_genome_id = ref_genome_tuple[0]
                    #all_inputs.append(os.path.join(species_folder, "processed" ,ref_genome_id, "mapped", f"{ind}_{ref_genome_id}_sorted.bam"))
                    #all_inputs.append(os.path.join(species_folder, "processed" ,ref_genome_id, "mapped", f"{ind}_{ref_genome_id}_sorted.bam.bai"))
                    #all_inputs.append(os.path.join(species_folder, "processed" ,ref_genome_id, "mapped", f"{ind}_{ref_genome_id}_sorted.rescaled.bam"))
                    #all_inputs.append(os.path.join(species_folder, "processed" ,ref_genome_id, "mapped", f"{ind}_{ref_genome_id}_sorted.rescaled.bam.bai"))

                    #all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "coverage", ind, f"{ind}_{ref_genome_id}_coverage_analysis.csv"))

                    #all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "coverage", f"{ref_genome_id}_combined_coverage_analysis.csv"))
                    #all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "coverage", f"{ref_genome_id}_combined_coverage_analysis_detailed.csv"))

                    #"{species}/results/{ref_genome}/endogenous/{individual}/{individual}_{ref_genome}.endogenous.csv"
                    #all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "endogenous", ind, f"{ind}_{ref_genome_id}.endogenous.csv"))

                    #"{species}/results/{ref_genome}/plots/endogenous_reads/{species}_endogenous_reads_pie_chart.png"
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "endogenous_reads", f"{species}_{ref_genome_id}_endogenous_reads_pie_chart.pdf"))

                    #depth plots
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "endogenous_reads", f"{species}_{ref_genome_id}_endogenous_reads_bar_chart.png"))
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_breadth_coverage_violin.png"))
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_depth_coverage_violin.png"))
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_depth_coverage_bar.png"))

                    # breadth plots
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_breadth_coverage_bins.png"))
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_depth_violin.png"))                
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_coverage_breadth_bar.png"))
                    all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_id, "plots", "coverage", f"{species}_{ref_genome_id}_individual_coverage_breadth_violin.png"))

                    all_inputs.append(os.path.join(species_folder, "processed" ,ref_genome_id, "consensus", f"{ind}_{ref_genome_id}", f"{ind}_{ref_genome_id}_consensus.fa.gz"))
                    #all_inputs.append(os.path.join(species_folder, "processed" ,ref_genome_id, "snp", f"{ind}_{ref_genome_id}", f"{ind}_{ref_genome_id}_snp.mafs.gz"))
            except Exception as e: 
                print(e)
                pass

        #for ref_genome_tuple in get_reference_genome_file_list_for_species(species):
            #{species}/results/{ref_genome}/endogenous/{ref_genome}_endogenous.csv
            #all_inputs.append(os.path.join(species_folder, "results" ,ref_genome_tuple[0], "endogenous", f"{ref_genome_tuple[0]}_endogenous.csv"))
        
    logging.info("Determined input for rule 'all':")
    for input in all_inputs:
        logging.info("\t" + "- " + input)

    return all_inputs
