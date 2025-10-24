import glob
import os

def expected_quality_filtered_files(wildcards):

    samples = get_sample_ids_for_species(wildcards.species)

    # find all samples that match the individual
    samples_of_individual = [f for f in samples if f.startswith(f"{wildcards.individual}")]
    
    if len(samples_of_individual) == 0:
        raise Exception(f"No raw read files found for individual {wildcards.individual}. Check that the individual ID is correct and that raw read files are present. Available samples: {samples}")
    
    # for each raw R1 file, generate the corresponding quality-filtered filename
    quality_filtered_files = []
    for sample in samples_of_individual:
        qf_file = f"{wildcards.species}/processed/reads/reads_quality_filtered/{sample}_quality_filtered.fastq.gz"
        quality_filtered_files.append(qf_file)
    
    return quality_filtered_files


# Rule: Merge quality-filtered reads by individual
rule merge_by_individual:
    input:
        expected_quality_filtered_files
    output:
        "{species}/processed/reads/reads_merged/{individual}.fastq.gz"
    message: 
        "Merging individual {wildcards.individual} of species {wildcards.species}."
    shell:
        """
        cat {input} > {output}
        """
