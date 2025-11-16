####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def standardize_reference_extension_to_fa_input_fasta(wc):

    # get ref paths for the given species
    reference_tuple = get_reference_file_list_for_species(wc.species)  # validate species and ref name

    reference_path = next((path for name, path in reference_tuple if name == wc.reference), None)

    if reference_path is None:
        raise ValueError(f"Reference {wc.reference} not found for species {wc.species}")

    return reference_path

####################################################
# Snakemake rules
####################################################

# Rule: Standardize reference extension to .fa
# 1) Normalize/standardize reference to .fa (symlink to avoid copying)
rule standardize_reference_extension_to_fa:
    input:
        fasta=standardize_reference_extension_to_fa_input_fasta
    output:
        fa="{species}/raw/ref/{reference}.fa"
    message:
        "Standardizing reference extension to .fa for {output.fa}"
    shell:
        "mv {input.fasta} {output.fa}"

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