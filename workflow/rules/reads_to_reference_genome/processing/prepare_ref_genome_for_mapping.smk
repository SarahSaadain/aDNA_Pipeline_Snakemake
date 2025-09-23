def find_reference_input(wc):
    """
    Return the first existing reference file for {species}/{ref_genome_name}
    among the allowed extensions.
    """

    # extensions we will accept for the original reference
    REF_EXTS = ["fa", "fasta", "fna", "fas", "fsa"]

    base = f"{wc.species}/raw/ref_genome/{wc.ref_genome_name}"
    for ext in REF_EXTS:
        p = f"{base}.{ext}"
        if os.path.exists(p):
            return p
    raise ValueError(
        f"No reference file found for {base}.{{" + ",".join(REF_EXTS) + "}}"
    )

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
    conda:
        "../../../envs/bwa.yaml"
    shell:
        "bwa index {input.fasta}"
