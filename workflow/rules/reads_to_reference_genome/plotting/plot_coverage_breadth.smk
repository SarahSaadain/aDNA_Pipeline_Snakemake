# Rule: Plot coverage breadth violin plot
rule coverage_violin_plot:
    input:
        "{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    output:
        report("{species}/results/{ref_genome}/plots/coverage/{species}_{ref_genome}_breadth_coverage_violin.png")
    message: "Plotting coverage breadth violin for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    params:
        species="{species}"
    log:
        "{species}/logs/{ref_genome}/plots/coverage/{species}_{ref_genome}_breadth_coverage_violin.log"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_coverage_breadth_violin.R"

# Rule: Plot coverage breadth bins
rule coverage_bins_plot:
    input:
        "{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    output:
        report("{species}/results/{ref_genome}/plots/coverage/{species}_{ref_genome}_breadth_coverage_bins.png")
    message: "Plotting coverage breadth bins for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    log:
        "{species}/logs/{ref_genome}/plots/coverage/{species}_{ref_genome}_breadth_coverage_bins.log"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_coverage_breadth_bins.R"

# Rule: Plot coverage breadth violin by individual
rule coverage_breadth_violin:
    input:
        "{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    output:
        report("{species}/results/{ref_genome}/plots/coverage/{species}_{ref_genome}_individual_coverage_breadth_violin.png")
    params:
        species="{species}"
    message: "Plotting coverage breadth violin for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    log:
        "{species}/logs/{ref_genome}/plots/coverage/{species}_{ref_genome}_individual_coverage_breadth_violin.log"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_coverage_breadth_by_individuals_violin.R"

# Rule: Plot coverage breadth bar by individual
rule coverage_breadth_bar:
    input:
        "{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    output:
        report("{species}/results/{ref_genome}/plots/coverage/{species}_{ref_genome}_individual_coverage_breadth_bar.png")
    message: "Plotting coverage breadth bar for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    params:
        species="{species}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_coverage_breadth_by_individuals_bar.R"
