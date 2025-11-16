#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(viridis)

# Access Snakemake variables
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]
species <- snakemake@params[["species"]]

# Read CSV
df <- read_csv(input_file)
if (nrow(df) == 0) stop("Input file is empty.")

# Extract individual
df <- df %>%
  mutate(individual = str_split(Filename, "_", simplify = TRUE)[,1])

# Aggregate average depth per individual (weighted mean)
df_summary <- df %>%
  group_by(individual) %>%
  summarise(avg_depth = weighted.mean(avg_depth, total_bases, na.rm = TRUE),
            .groups = "drop")

bar_plot <- ggplot(df_summary, aes(x = reorder(individual, -avg_depth),
                                   y = avg_depth, fill = individual)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  ylab("Avg. Depth") +
  xlab("Individual") +
  ggtitle(paste0("Average Depth per Individual (Aggregated): ", species)) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  scale_fill_viridis_d()

ggsave(output_file, plot = bar_plot, width = 12, height = 6, dpi = 300)
