####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def get_final_bam_input_bam(wildcards):

    species = wildcards["species"]
    reference = wildcards["reference"]
    individual = wildcards["individual"]

    # if damage analysis is enabled, return the rescaled bam
    # if deduplication is enabled, return the dedupped bam
    # if map_reads_to_reference is enabled, return the sorted bam

    if config.get("pipeline", {}).get("reference_processing", {}).get("damage_rescaling", {}).get("execute", True) == True:
        return f"{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped_rescaled.bam"

    if config.get("pipeline", {}).get("reference_processing", {}).get("deduplication", {}).get("execute", True) == True:
        return f"{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped.bam"

    return f"{species}/processed/{reference}/mapped/{individual}_{reference}_sorted.bam"

def get_final_bam_input_bai(wildcards):

    return f"{get_final_bam_input_bam(wildcards)}.bai"

####################################################
# Snakemake rules
####################################################

rule get_final_bam:
    input:
        bam = get_final_bam_input_bam,
        bai = get_final_bam_input_bai
    output:
        bam = "{species}/processed/{reference}/mapped/{individual}_{reference}_final.bam",
        bai = "{species}/processed/{reference}/mapped/{individual}_{reference}_final.bam.bai"
    
    message:
        "Getting final bam for {wildcards.individual} of {wildcards.species}."
    shell:
        # The final bam is just a copy of the bam that is output from the last step of the reference 
        # processing pipeline (either damage rescaling, deduplication, or mapping)
        # We use copy instead of move. If we would use move here, other rules will not know that the input bam 
        # will be renamed. These other rules will then fail because the input bam is not found (as it has been renamed). 
        # By using copy, we ensure that the input bam is still available for other rules that need it. 
        # The input bam will be deleted at the end of the pipeline when the intermediate files are cleaned up.
        """
        echo "Copying final bam and bai for {wildcards.individual} mapped to {wildcards.reference}..."
        echo "Input bam: {input.bam}"
        echo "Output bam: {output.bam}"
        cp {input.bam} {output.bam}
        cp {input.bai} {output.bai}
        echo "Done copying final bam and bai for {wildcards.individual} mapped to {wildcards.reference}."
        """