#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(viridis)

species <- snakemake@params[["species"]]
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

df <- read_csv(input_file)
if (nrow(df) == 0) stop("Input file is empty.")

df <- df %>%
  mutate(individual = str_split(Filename, "_", simplify = TRUE)[,1])

# Bar plot function: weighted overall percent covered
df_summary <- df %>%
  group_by(individual) %>%
  summarise(
    total_covered = sum(covered_bases, na.rm = TRUE),
    total_bases = sum(total_bases, na.rm = TRUE),
    overall_percent_covered = ifelse(total_bases > 0, 100 * total_covered / total_bases, 0),
    .groups = "drop"
  )

bar_plot <- ggplot(df_summary, aes(x = reorder(individual, -overall_percent_covered),
                                   y = overall_percent_covered, fill = individual)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Overall Percent Covered") +
  xlab("Individual") +
  ggtitle(paste0("Overall Breadth of Coverage per Individual: ", species)) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  scale_fill_viridis_d()

ggsave(output_file, plot = bar_plot, width = 12, height = 6, dpi = 300)
