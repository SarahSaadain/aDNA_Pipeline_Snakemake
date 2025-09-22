rule plot_endogenous_reads_pie:
    input:
        "{species}/results/{ref_genome}/endogenous/{ref_genome}_endogenous.csv"
    output:
        report("{species}/results/{ref_genome}/plots/endogenous_reads/{species}_{ref_genome}_endogenous_reads_pie_chart.pdf")
    message: "Plotting endogenous reads pie chart for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_endogenous_reads_pie.R"

rule plot_endogenous_reads_bar:
    input:
        "{species}/results/{ref_genome}/endogenous/{ref_genome}_endogenous.csv"
    output:
        report("{species}/results/{ref_genome}/plots/endogenous_reads/{species}_{ref_genome}_endogenous_reads_bar_chart.png")
    message: "Plotting endogenous reads bar chart for species {wildcards.species} and reference genome {wildcards.ref_genome}"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference_genome/plotting/plot_endogenous_reads_bar.R"