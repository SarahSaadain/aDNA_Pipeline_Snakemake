library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)

plot_endogenous_reads_bar <- function(source_file, output_file) {
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Read the TSV file with pipeline stages
  df <- read.table(source_file, sep = "\t", header = TRUE) %>%
    mutate(
      individual = sub("_.*", "", individual),
      percent_endogenous = (mapped_endogenous_reads / raw_reads) * 100,
      percent_after_dup_removal = (endogenous_duplicates_removed / raw_reads) * 100,
      percent_duplicates_removed = percent_endogenous - percent_after_dup_removal
    )
  
  # Create long format for stacked bars showing duplicates
  df_long <- df %>%
    select(individual, percent_after_dup_removal, percent_duplicates_removed) %>%
    pivot_longer(
      cols = -individual,
      names_to = "category",
      values_to = "percent"
    ) %>%
    mutate(
      category = factor(category,
                        levels = c("percent_after_dup_removal",
                                  "percent_duplicates_removed"))
    )
  
  p <- ggplot(df_long, aes(x = individual, y = percent, fill = category)) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    theme_bw() +
    labs(x = "Individual", y = "Percentage of Endogenous Reads", title = "Endogenous Reads including Duplicate Removal") + 
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black")) +
    scale_fill_manual(
      name = "Category",
      values = c(
        "percent_after_dup_removal" = "#2E86AB",
        "percent_duplicates_removed" = "#A23B72"
      ),
      labels = c(
        "percent_after_dup_removal" = "Final Endogenous",
        "percent_duplicates_removed" = "Duplicates Removed"
      )
    )
  
  
  ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
}

if (exists("snakemake")) {
  source_file <- snakemake@input[[1]]
  output_file <- snakemake@output[[1]]
  plot_endogenous_reads_bar(source_file, output_file)
}