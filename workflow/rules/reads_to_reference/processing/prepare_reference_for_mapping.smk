####################################################
# Snakemake rules
####################################################
rule standardize_reference_extension_to_fa:
    output:
        fa="{species}/raw/ref/{reference}.fa"
    message:
        "Ensuring reference {wildcards.reference} for {wildcards.species} is standardized to .fa"
    conda:
        "../../../envs/python.yaml",
    run:
        # we use python to rename the file if necessary
        import os
        import shutil

        # Get the list of reference tuples (sanitized_name, full_path)
        # the full_path contains the original file path
        reference_tuples = get_reference_file_list_for_species(wildcards.species)

        # Find the path corresponding to the sanitized reference name
        ref_path = next((path for name, path in reference_tuples if name == wildcards.reference), None)

        if ref_path is None:
            logger.error(f"Reference {wildcards.reference} not found for species {wildcards.species}")
            logger.error(f"Available references: {reference_tuples}")
            raise ValueError(f"Reference {wildcards.reference} not found for species {wildcards.species}")

        if not os.path.exists(ref_path):
            logger.error(f"Reference file {ref_path} does not exist.")
            raise FileNotFoundError(f"Reference file {ref_path} does not exist.")

        # Create the output folder if it doesn't exist
        os.makedirs(os.path.dirname(output.fa), exist_ok=True)

        # Only rename if the standardized file doesn't already exist
        if not os.path.exists(output.fa):
            # Use symlink if you don't want to copy the file
            os.rename(ref_path, output.fa)
            logger.info(f"Reference {ref_path} renamed to {output.fa}")
        else:
            logger.info(f"Reference {output.fa} already exists, skipping.")


# Rule: Index reference with BWA
# 2) BWA index on the standardized .fa
rule index_reference_for_mapping:
    input:
        "{species}/raw/ref/{reference}.fa"
    output:
        multiext("{species}/raw/ref/{reference}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    message: "Indexing reference {wildcards.reference} with BWA"
    log:
        "{species}/logs/{reference}/index/{reference}_bwa_index.log"
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v7.2.0/bio/bwa/index"

rule index_reference_with_samtools:
    input:
       "{species}/raw/ref/{reference}.fa"
    output:
        "{species}/raw/ref/{reference}.fa.fai"
    log:
        "{species}/raw/ref/{reference}.fa.fai.log"
    params:
        extra="",
    wrapper:
        "v7.6.0/bio/samtools/faidx"