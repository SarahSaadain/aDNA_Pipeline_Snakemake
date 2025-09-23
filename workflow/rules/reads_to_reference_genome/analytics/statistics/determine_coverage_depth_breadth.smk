# Rule: Calculate coverage depth using samtools
rule samtools_depth:
    input:
        bams=["{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam"]
    output:
        "{species}/processed/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_depth.tsv"
    log:
        "{species}/processed/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_depth.log"
    message: "Calculating coverage depth for {input.bams}"
    params:
        # optional bed file passed to -b
        extra="-aa",  # optional additional parameters as string
    wrapper:
        "v7.2.0/bio/samtools/depth"

# Rule: Analyze coverage depth and breadth
rule samtools_coverage_analysis:
    input:
        depth_txt="{species}/processed/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_depth.tsv"
    output:
        analysis="{species}/results/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_coverage_analysis.csv"
    message: "Analyzing coverage depth and breadth for {input.depth_txt}"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/analyze_samtools_depth_file.py"

# Rule: Combine coverage analysis files
rule combine_coverage:
    input:
        analysis = lambda wildcards: expand(
            "{species}/results/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_coverage_analysis.csv",
            species=wildcards.species,
            ref_genome=wildcards.ref_genome,
            individual=get_individuals_for_species(wildcards.species),
        )
    output:
        combined="{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis.csv",
        detailed="{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    message: "Combining coverage analysis files for species {params.species}"
    params:
        species="{species}"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/combine_analyzed_depth_breadth_files.py"
