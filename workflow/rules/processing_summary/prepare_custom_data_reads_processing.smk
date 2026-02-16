####################################################
# Snakemake rules
####################################################

def prepare_custom_data_reads_processing_dedup(wildcards):

    if config.get("pipeline", {}).get("reference_processing", {}).get("execute", True) == False:
        return None
    
    if config.get("pipeline", {}).get("reference_processing", {}).get("deduplication", {}).get("execute", True) == True:
        return f"{wildcards.species}/results/{wildcards.reference}/analytics/{wildcards.individual}/dedup/{wildcards.individual}_{wildcards.reference}_final.dedup.json"
    else:
        return None

def prepare_custom_data_reads_processing_endogenous(wildcards):
    
    if config.get("pipeline", {}).get("reference_processing", {}).get("execute", True) == True:
        #"{species}/results/{reference}/analytics/{individual}/endogenous/{individual}_{reference}.endogenous.csv"
        return f"{wildcards.species}/results/{wildcards.reference}/analytics/{wildcards.individual}/endogenous/{wildcards.individual}_{wildcards.reference}.endogenous.csv"
    else:
        return None

####################################################
# Snakemake rules
####################################################

rule prepare_custom_data_reads_processing:
    input:
        reads = "{species}/results/reads/statistics/{species}_reads_counts.csv",
        endogenous = prepare_custom_data_reads_processing_endogenous,
        dedup = prepare_custom_data_reads_processing_dedup
    output:
        "{species}/results/summary/{individual}/multiqc_custom_content/{individual}_{reference}_reads_processing_summary.tsv",
    params:
        individual="{individual}",
        reference="{reference}",
    script:
        "../../scripts/processing_summary/prepare_custom_data_reads_processing.py"
        