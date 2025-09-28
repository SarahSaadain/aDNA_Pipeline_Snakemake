import glob
import os

def expected_quality_filtered_files(wildcards):
    # folder with raw reads
    raw_folder = f"{wildcards.species}/raw/reads"
    
    # find all raw R1 files for this individual
    raw_files_r1 = glob.glob(os.path.join(raw_folder, f"{wildcards.individual}_*R1*.fastq.gz"))
    
    if len(raw_files_r1) == 0:
        raise Exception(f"No raw R1 files found for individual {wildcards.individual}")
    
    # for each raw R1 file, generate the corresponding quality-filtered filename
    quality_filtered_files = []
    for raw_file in raw_files_r1:
        filename = os.path.basename(raw_file)
        # get the part between individual and _R1
        rest = filename[len(wildcards.individual)+1:filename.index("_R1")]
        qf_file = f"{wildcards.species}/processed/reads/reads_quality_filtered/{wildcards.individual}_{rest}_quality_filtered.fastq.gz"
        quality_filtered_files.append(qf_file)
    
    return quality_filtered_files


# Rule: Merge quality-filtered reads by individual
rule merge_by_individual:
    input:
        expected_quality_filtered_files
    output:
        temp("{species}/processed/reads/reads_merged/{individual}.fastq.gz")
    message: 
        "Merging individual {wildcards.individual} of species {wildcards.species}."
    shell:
        """
        cat {input} > {output}
        """
