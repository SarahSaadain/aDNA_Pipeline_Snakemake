rule plot_depth_violin:
    input:
        depth_file = "{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    output:
        plot = report("{species}/results/{ref_genome}/plots/coverage/{species}_{ref_genome}_depth_violin.png")
    message: "Plotting depth violin for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    params:
        species="{species}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_coverage_depth_violin.R"

rule depth_coverage_violin:
    input:
        "{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    output:
        report("{species}/results/{ref_genome}/plots/coverage/{species}_{ref_genome}_individual_depth_coverage_violin.png")
    message: "Plotting depth coverage violin for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    params:
        species="{species}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_coverage_depth_by_individuals_violin.R"

rule depth_coverage_bar:
    input:
        "{species}/results/{ref_genome}/coverage/{ref_genome}_combined_coverage_analysis_detailed.csv"
    output:
        report("{species}/results/{ref_genome}/plots/coverage/{species}_{ref_genome}_individual_depth_coverage_bar.png")
    message: "Plotting depth coverage bar for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    params:
        species="{species}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_coverage_depth_by_individuals_bar.R"
