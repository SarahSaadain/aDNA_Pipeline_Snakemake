####################################################
# Python helper functions for rules
# Naming of functions: <rule_name>_<rule_parameter>[_<rule_subparameter>]>
####################################################

def helper_get_bam_for_damage_analysis(wildcards):

    species = wildcards["species"]
    reference_id = wildcards["reference"]
    ind = wildcards["individual"]

    # if deduplication is enabled, return the dedupped bam
    # if map_reads_to_reference is enabled, return the sorted bam

    if config.get("pipeline", {}).get("reference_processing", {}).get("deduplication", {}).get("execute", True) == True:
        return os.path.join(species, "processed" ,reference_id, "mapped", f"{ind}_{reference_id}_sorted_dedupped.bam")

    return os.path.join(species, "processed" ,reference_id, "mapped", f"{ind}_{reference_id}_sorted.bam")


def analyze_damageprofile_input_bam(wildcards):
    return helper_get_bam_for_damage_analysis(wildcards)
    
def analyze_mapdamage_and_rescale_bam_input_bam(wildcards):
    return helper_get_bam_for_damage_analysis(wildcards)

def analyze_mapdamage_and_rescale_bam_input_bam_index(wildcards):
    return f"{analyze_damageprofile_input_bam(wildcards)}.bai"

#"{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped.bam",

####################################################
# Snakemake rules
####################################################

rule analyze_damageprofile:
    input:
        bam = analyze_damageprofile_input_bam,
        ref = "{species}/raw/ref/{reference}.fa",
        ref_index = "{species}/raw/ref/{reference}.fa.fai"
    output:
        damage_profile = directory("{species}/results/{reference}/damage/damageprofile/{individual}")
    message:
        "Generate damage profile for {wildcards.individual} mapped to {wildcards.reference}",
    conda:
        "../../../envs/damage_profiler.yaml"
    shell:
        """
        mkdir -p {output.damage_profile}
        damageprofiler -i {input.bam} -r {input.ref} -o {output.damage_profile}
        """

# Rule: Analyze DNA damage and rescale BAM using mapDamage2
rule analyze_mapdamage_and_rescale_bam:
    input:
        bam = analyze_mapdamage_and_rescale_bam_input_bam,
        bam_index = analyze_mapdamage_and_rescale_bam_input_bam_index,
        ref = "{species}/raw/ref/{reference}.fa"
    output:
        log = "{species}/results/{reference}/damage/mapdamage/{individual}/Runtime_log.txt",
        GtoA3p = "{species}/results/{reference}/damage/mapdamage/{individual}/3pGtoA_freq.txt",
        CtoT5p = "{species}/results/{reference}/damage/mapdamage/{individual}/5pCtoT_freq.txt",
        dnacomp = "{species}/results/{reference}/damage/mapdamage/{individual}/dnacomp.txt",
        frag_misincorp = "{species}/results/{reference}/damage/mapdamage/{individual}/Fragmisincorporation_plot.pdf",
        #len = "{species}/results/{reference}/damage/{individual}/Length_plot.pdf",
        lg_dist = "{species}/results/{reference}/damage/mapdamage/{individual}/lgdistribution.txt",
        misincorp = "{species}/results/{reference}/damage/mapdamage/{individual}/misincorporation.txt",
        rescaled_bam = temp("{species}/results/{reference}/damage/mapdamage/{individual}/{individual}_{reference}.bam")
    params:
        extra="--merge-reference-sequences --rescale",  # optional parameters for mapdamage2 (except -i, -r, -d, --rescale)
    message:
        "Analyze damage and rescale {input.bam}",
    benchmark:
        "benchmark/{species}_{individual}_{reference}_mapdamage2.benchmark.json",
    log:
        "{species}/logs/{reference}/damage/mapdamage/{individual}/mapdamage2.log",
    wrapper:
        "v7.2.0/bio/mapdamage2"

# Rule: Sort rescaled BAM file
# 3 Sort BAM
rule sort_rescaled_bam:
    input:
        "{species}/results/{reference}/damage/mapdamage/{individual}/{individual}_{reference}.bam"
    output:
        "{species}/results/{reference}/damage/mapdamage/{individual}/{individual}_{reference}_sorted.bam"
    log:
        "{species}/logs/{reference}/damage/mapdamage/{individual}/{individual}_{reference}_sorted.bam.log"
    message:
        "Sort rescaled BAM for {input}",
    threads: 10
    wrapper:
        "v7.5.0/bio/samtools/sort"

# Rule: Index rescaled BAM file
# 4 Index BAM
rule index_rescaled_bam:
    input:
        "{species}/results/{reference}/damage/mapdamage/{individual}/{individual}_{reference}_sorted.bam"
    output:
        "{species}/results/{reference}/damage/mapdamage/{individual}/{individual}_{reference}_sorted.bam.bai"
    log:
        "{species}/logs/{reference}/damage/mapdamage/{individual}/{individual}_{reference}_sorted.bam.bai.log"
    message:
        "Index rescaled BAM for {input}",
    threads: 10
    wrapper:
        "v7.5.0/bio/samtools/index"

# Rule: Move rescaled BAM and index to processed directory
rule move_rescaled_bam:
    input:
        sorted_bam="{species}/results/{reference}/damage/mapdamage/{individual}/{individual}_{reference}_sorted.bam",
        bam_index ="{species}/results/{reference}/damage/mapdamage/{individual}/{individual}_{reference}_sorted.bam.bai"
    output:
        sorted_bam="{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped_rescaled.bam",
        bam_index ="{species}/processed/{reference}/mapped/{individual}_{reference}_sorted_dedupped_rescaled.bam.bai"
    message:
        "Move rescaled BAM and index to processed directory for {input.sorted_bam}",
    shell:
        """
        mv {input.sorted_bam} {output.sorted_bam} 
        mv {input.bam_index} {output.bam_index} 
        """

