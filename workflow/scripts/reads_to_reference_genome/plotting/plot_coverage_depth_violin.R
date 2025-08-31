#!/usr/bin/env Rscript

# Load libraries
library(ggplot2)
library(dplyr)
library(scales)

# Function to generate depth coverage violin plot
plot_depth_coverage_violin <- function(df_depth, species) {
  df_depth$species <- species

  p <- ggplot(df_depth, aes(x = factor(species), y = avg_depth)) +
    scale_y_continuous(
      trans = "log10",
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = comma
    ) +
    geom_violin(scale = "width", fill = "#209557") +
    theme_bw() +
    ylab("Avg. Depth") +
    xlab("Species") +
    ggtitle("Distribution of Average Depth") +
    theme(
      axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.position = "none"
    )
  return(p)
}

# Save plot function
save_plot <- function(plot, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(output_file, plot = plot, width = 10, height = 6)
  cat("✅ Saved:", output_file, "\n")
}

# === Main ===

# Snakemake mode
species <- snakemake@params[["species"]]
depth_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]


cat("🎯 Plotting depth coverage violin for species:", species, "\n")
cat("Input file:", depth_file, "\n")
cat("Output file:", output_file, "\n")

# Read depth file
df_depth <- read.csv(depth_file, header = TRUE)

if (nrow(df_depth) == 0) {
  stop("❌ Input depth file is empty.")
}

# Generate and save plot
plot_violin <- plot_depth_coverage_violin(df_depth, species)
save_plot(plot_violin, output_file)
