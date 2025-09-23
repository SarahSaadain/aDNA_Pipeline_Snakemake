# Rule: Determine endogenous reads from BAM stats
rule endogenous_reads:
    input:
        stats="{species}/results/{ref_genome}/statistics/{individual}/{individual}_{ref_genome}.bam.stats"
    output:
        csv= "{species}/results/{ref_genome}/endogenous/{individual}/{individual}_{ref_genome}.endogenous.csv"
    message: "Determining endogenous reads for {input.stats}"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/parse_endogenous_from_stats.py"

# Rule: Combine endogenous reads for all individuals
rule combine_endogenous_reads:
    input:
        lambda wildcards: expand(
            "{species}/results/{ref_genome}/endogenous/{individual}/{individual}_{ref_genome}.endogenous.csv",
            species=wildcards.species,
            ref_genome=wildcards.ref_genome,
            individual=get_individuals_for_species(wildcards.species),
        )
    output:
        "{species}/results/{ref_genome}/endogenous/{ref_genome}_endogenous.csv"
    message: "Combining endogenous reads for species {wildcards.species}"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/combine_endogenous_reads.py"
