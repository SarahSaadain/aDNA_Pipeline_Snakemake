rule plot_read_counts:
    input:
        "{species}/results/qualitycontrol/statistics/{species}_reads_processing.csv"
    output:
        report("{species}/results/plots/{species}_read_count_comparison.png")
    message: "Plotting read counts comparison for species {wildcards.species}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/raw_reads/plotting/plot_comparison_reads_before_after_processing.R"

rule plot_read_counts_comparison_by_individual:
    input:
        "{species}/results/qualitycontrol/statistics/{species}_reads_processing.csv"
    output:
         report("{species}/results/plots/{species}_read_count_comparison_by_individual.png")
    message: "Plotting read counts comparison per individual for species {wildcards.species}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/raw_reads/plotting/plot_comparison_reads_before_after_processing_per_individual.R"
