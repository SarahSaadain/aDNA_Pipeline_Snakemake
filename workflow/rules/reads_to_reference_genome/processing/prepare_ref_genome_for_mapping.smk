def find_reference_input(wc):

    # get ref genome paths for the given species
    ref_genomes_tuple = get_reference_genome_file_list_for_species(wc.species)  # validate species and ref genome name

    ref_genome_path = next((path for name, path in ref_genomes_tuple if name == wc.ref_genome_name), None)

    if ref_genome_path is None:
        raise ValueError(f"Reference genome {wc.ref_genome_name} not found for species {wc.species}")

    return ref_genome_path

# Rule: Standardize reference genome extension to .fa
# 1) Normalize/standardize reference to .fa (symlink to avoid copying)
rule standardize_reference_extension_to_fa:
    input:
        fasta=find_reference_input
    output:
        fa="{species}/raw/ref_genome/{ref_genome_name}.fa"
    message:
        "Standardizing reference genome extension to .fa for {output.fa}"
    shell:
        "mv {input.fasta} {output.fa}"

# Rule: Index reference genome with BWA
# 2) BWA index on the standardized .fa
rule bwa_index:
    input:
        fasta="{species}/raw/ref_genome/{ref_genome_name}.fa"
    output:
        amb="{species}/raw/ref_genome/{ref_genome_name}.fa.amb",
        ann="{species}/raw/ref_genome/{ref_genome_name}.fa.ann",
        bwt="{species}/raw/ref_genome/{ref_genome_name}.fa.bwt",
        pac="{species}/raw/ref_genome/{ref_genome_name}.fa.pac",
        sa="{species}/raw/ref_genome/{ref_genome_name}.fa.sa"
    message: "Indexing reference genome {input.fasta} with BWA"
    log:
        "{species}/logs/ref_genome_name/index/{ref_genome_name}_bwa_index.log"
    conda:
        "../../../envs/bwa.yaml"
    shell:
        "bwa index {input.fasta} > {log} 2>&1"
