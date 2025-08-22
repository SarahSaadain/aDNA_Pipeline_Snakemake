library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tools)

# Violin plot: per-scaffold avg depth distribution by individual
plot_depth_coverage_violin <- function(df) {
  ggplot(df, aes(x = factor(individual), y = avg_depth, fill = individual)) +
    scale_y_continuous(
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::comma
    ) +
    geom_violin(scale = "width", trim = TRUE) +
    theme_bw() +
    ylab("Avg. Depth (log10)") +
    xlab("Individual") +
    ggtitle("Per-Scaffold Average Depth per Individual") +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black")) +
    scale_fill_viridis_d()
}

# Bar plot: average depth per individual (aggregated)
plot_depth_coverage_bar <- function(df) {
  # Aggregate average depth per individual (simple mean or weighted mean)
  df_summary <- df %>%
    group_by(individual) %>%
    summarise(avg_depth = weighted.mean(avg_depth, total_bases, na.rm = TRUE)) %>%
    ungroup()

  ggplot(df_summary, aes(x = reorder(individual, -avg_depth), y = avg_depth, fill = individual)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ylab("Avg. Depth") +
    xlab("Individual") +
    ggtitle("Average Depth per Individual (Aggregated)") +
    theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
    scale_fill_viridis_d()
}

# Command-line args: <species> <input_file.csv> <output_folder>
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("❌ Usage: Rscript script.R <species> <input_file.csv> <output_folder>")
}

species <- args[1]
input_file <- args[2]
output_folder <- args[3]

cat("🎯 Plotting coverage depth for ", species, "\n")
cat("Input file: ", input_file, "\n")
cat("Output folder: ", output_folder, "\n")

# Read detailed per-scaffold CSV
df <- read_csv(input_file)

if (nrow(df) == 0) {
  stop("Input file is empty.")
}

# Add individual column from Filename (e.g., 'Phor01.fastq_marker_...' -> 'Phor01')
df <- df %>%
  mutate(individual = str_split(Filename, "_", simplify = TRUE)[,1]) %>%
  rename(avg_depth = avg_depth, total_bases = total_bases)  # in case column names differ

# Ensure output directory exists
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
  cat("Created output directory:", output_folder, "\n")
}

# Create and save plots
violin_plot <- plot_depth_coverage_violin(df)
bar_plot <- plot_depth_coverage_bar(df)

violin_file <- file.path(output_folder, paste0(species, "_violin_per_scaffold_avg_depth.png"))
bar_file <- file.path(output_folder, paste0(species, "_barplot_aggregated_avg_depth.png"))

ggsave(violin_file, plot = violin_plot, width = 12, height = 6)
ggsave(bar_file, plot = bar_plot, width = 12, height = 6)

cat("✅ Saved violin plot to:", violin_file, "\n")
cat("✅ Saved bar plot to:", bar_file, "\n")