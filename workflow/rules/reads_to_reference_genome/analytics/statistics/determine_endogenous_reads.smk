import scripts.setup_snakemake as setup
import input_manager as im


rule endogenous_reads:
    input:
        stats="{species}/results/{ref_genome}/statistics/{individual}/{individual}_{ref_genome}.bam.stats"
    output:
        csv= "{species}/results/{ref_genome}/endogenous/{individual}/{individual}_{ref_genome}.endogenous.csv"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/parse_endogenous_from_stats.py"

rule combine_endogenous_reads:
    input:
        lambda wildcards: expand(
            "{species}/results/{ref_genome}/endogenous/{individual}/{individual}_{ref_genome}.endogenous.csv",
            species=wildcards.species,
            ref_genome=wildcards.ref_genome,
            individual=im.get_individuals_for_species(wildcards.species),
        )
    output:
        "{species}/results/{ref_genome}/endogenous/{ref_genome}_endogenous.csv"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/combine_endogenous_reads.py"
