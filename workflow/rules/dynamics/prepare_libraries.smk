####################################################
# Snakemake rules
####################################################
rule prepare_feature_library:
    input:
        # Get the TE library for the species from the config file
        "{species}/raw/dynamics/feature_library/{feature_library}.fasta"
    output:
        temp("{species}/processed/dynamics/lib/feature_library/{feature_library}.suffixed.fasta")
    message: "Preparing TE library for {wildcards.species}"
    shell:
        # remove trailing whitespace from headers and append _fle to each header
        """
        sed -E '/^>/ s/[[:space:]]//g; /^>/ s/(_fle)?$/_fle/' {input} > {output}
        """

rule prepare_scg_library:
    input:
        "{species}/raw/dynamics/scg/{scg_library}.fasta"
    output:
        temp("{species}/processed/dynamics/lib/scg/{scg_library}.suffixed.fasta")
    message: "Preparing SCG library for {wildcards.species}"
    shell:
        # remove trailing whitespace from headers and append _scg to each header
        """
        sed -E '/^>/ s/[[:space:]]//g; /^>/ s/(_scg)?$/_scg/' {input} > {output}
        """
        
rule combine_scg_and_ref_library:
    input:
        scg= lambda wildcards: f"{wildcards.species}/processed/dynamics/lib/scg/{get_scg_library_ids_for_species(wildcards.species)[0]}.suffixed.fasta",
        fl="{species}/processed/dynamics/lib/feature_library/{feature_library}.suffixed.fasta"
    output:
        library="{species}/processed/dynamics/lib/{feature_library}_and_scg.suffixed.fasta"
    message: "Concatenating SCG and Feature libraries for {wildcards.species}"
    shell:
        """
        cat {input.scg} {input.fl} > {output.library}
        """