####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def get_final_bam_input_bam(wildcards):

    species = wildcards["species"]
    reference_id = wildcards["reference"]
    ind = wildcards["individual"]

    # if damage analysis is enabled, return the rescaled bam
    # if deduplication is enabled, return the dedupped bam
    # if map_reads_to_reference is enabled, return the sorted bam

    if config.get("pipeline", {}).get("reference_processing", {}).get("damage_rescaling", {}).get("execute", True) == True:
        return os.path.join(species, "processed" ,reference_id, "mapped", f"{ind}_{reference_id}_sorted_dedupped_rescaled.bam")
    
    if config.get("pipeline", {}).get("reference_processing", {}).get("deduplication", {}).get("execute", True) == True:
        return os.path.join(species, "processed" ,reference_id, "mapped", f"{ind}_{reference_id}_sorted_dedupped.bam")

    return os.path.join(species, "processed" ,reference_id, "mapped", f"{ind}_{reference_id}_sorted.bam")

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
    shell:
        """
        mv {input.bam} {output.bam}
        mv {input.bai} {output.bai}
        """