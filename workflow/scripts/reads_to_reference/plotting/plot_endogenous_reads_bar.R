library(dplyr)
library(ggplot2)
library(viridis)

plot_endogenous_reads_bar <- function(source_file, output_file) {
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Read the combined CSV file
  df <- read.table(source_file, sep =",", header = TRUE) %>%
    mutate(
      non_endogenous = total_reads - mapped_reads,
      percent_endogenous = proportion * 100
    )
  
  # The 'Individual' column is now loaded directly from the input file, 
  # so no need to extract it from 'Filename'.
  
  p <- ggplot(df, aes(x = individual, y = percent_endogenous, fill = individual)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    # Updated X-axis label to reflect the use of the 'Individual' column
    labs(x = "Individual", y = "Percentage of Endogenous Reads", title = "Endogenous Reads") + 
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black")) +
    scale_fill_viridis_d()
  
  ggsave(output_file, plot = p, width = 6, height = 6, dpi = 300)
}

if (exists("snakemake")) {
  source_file <- snakemake@input[[1]]
  output_file <- snakemake@output[[1]]
  plot_endogenous_reads_bar(source_file, output_file)
}