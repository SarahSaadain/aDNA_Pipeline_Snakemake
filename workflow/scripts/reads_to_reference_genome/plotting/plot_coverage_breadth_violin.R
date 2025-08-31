library(ggplot2)
library(dplyr)
library(readr)

plot_coverage_breadth_violin <- function(df, species, output_file) {
  
  df$species <- species
  
  plot_breadth <- ggplot(df, aes(x = factor(species), y = percent_covered)) +
    geom_violin(scale = "width") +
    theme_bw() +
    ylab("Percent Covered") +
    xlab("Species") +
    ggtitle("Distribution of Percent Covered") +
    theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          legend.position = "none")
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(output_file, plot = plot_breadth, width = 10, height = 6)
  cat("✅ Saved:", output_file, "\n")
}

# Snakemake integration
if (exists("snakemake")) {
  species <- snakemake@params[["species"]]
  input_file <- snakemake@input[[1]]
  output_file <- snakemake@output[[1]]
  
  df <- read.csv(input_file, header = TRUE)
  plot_coverage_breadth_violin(df, species, output_file)
}
