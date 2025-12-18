####################################################
# Snakemake rules
####################################################

# Rule: Plot endogenous reads bar chart
rule plot_mapped_reads_endogenous_bar:
    input:
        "{species}/results/{reference}/analytics/{species}/endogenous/{reference}_endogenous.csv"
    output:
        plot = "{species}/results/{reference}/plots/endogenous_reads/{species}_{reference}_endogenous_reads_bar_chart.png"
    message: "Plotting endogenous reads bar chart for species {wildcards.species} and reference {wildcards.reference}"
    log:
        "{species}/logs/{reference}/plots/endogenous_reads/{species}_{reference}_endogenous_reads_bar_chart.log"
    conda:
        "../../../envs/r_plot.yaml"
    script:
        "../../../scripts/reads_to_reference/plotting/plot_endogenous_reads_bar.R"