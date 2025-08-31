import scripts.setup_snakemake as setup
import input_manager as im

rule samtools_depth:
    input:
        bams=["{species}/processed/{ref_genome}/mapped/{individual}_{ref_genome}_sorted.bam"]
    output:
        "{species}/processed/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_depth.tsv"
    log:
        "{species}/processed/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_depth.log"
    params:
        # optional bed file passed to -b
        extra="",  # optional additional parameters as string
    wrapper:
        "v7.2.0/bio/samtools/depth"

rule samtools_coverage_analysis:
    input:
        depth_txt="{species}/processed/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_depth.tsv"
    output:
        analysis="{species}/results/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_coverage_analysis.csv"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/analyze_samtools_depth_file.py"

rule combine_coverage:
    input:
        analysis = lambda wildcards: expand(
            "{species}/results/{ref_genome}/coverage/{individual}/{individual}_{ref_genome}_coverage_analysis.csv",
            species=wildcards.species,
            ref_genome=wildcards.ref_genome,
            individual=im.get_individuals_for_species(wildcards.species),
        )
    output:
        combined="{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis.csv",
        detailed="{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    params:
        species="{species}"
    script:
        "../../../../scripts/reads_to_reference_genome/analytics/statistics/combine_analyzed_depth_breadth_files.py"
