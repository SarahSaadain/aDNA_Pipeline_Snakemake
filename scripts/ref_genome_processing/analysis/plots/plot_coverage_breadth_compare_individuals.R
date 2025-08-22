library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tools)

# Violin plot function: percent covered per scaffold
plot_breadth_violin <- function(df) {
  ggplot(df, aes(x = factor(individual), y = percent_covered, fill = individual)) +
    geom_violin(scale = "width", trim = TRUE) +
    theme_bw() +
    ylab("Percent Covered (per scaffold)") +
    xlab("Individual") +
    ggtitle("Distribution of Coverage Breadth per Scaffold per Individual") +
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

# Bar plot function: overall percent covered per individual
plot_breadth_bar <- function(df) {
  # Weighted breadth calculation: total covered / total bases
  df_summary <- df %>%
    group_by(individual) %>%
    summarise(
      total_covered = sum(covered_bases, na.rm = TRUE),
      total_bases = sum(total_bases, na.rm = TRUE),
      overall_percent_covered = ifelse(total_bases > 0, 100 * total_covered / total_bases, 0)
    ) %>%
    ungroup()

  ggplot(df_summary, aes(x = reorder(individual, -overall_percent_covered), y = overall_percent_covered, fill = individual)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    ylab("Overall Percent Covered") +
    xlab("Individual") +
    ggtitle("Overall Breadth of Coverage per Individual") +
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
  stop("❌ Usage: Rscript breadth_plot_from_detailed.R <species> <input_file.csv> <output_folder>")
}

species <- args[1]
input_file <- args[2]
output_folder <- args[3]

cat("🎯 Plotting coverage breadth for species ", species)
cat("Input file: ", input_file)
cat("Target folder: ", output_folder)

# Load data
df <- read_csv(input_file)

if (nrow(df) == 0) {
  stop("Input file is empty.")
}

# Extract individual from filename
df <- df %>%
  mutate(individual = str_split(Filename, "_", simplify = TRUE)[,1]) %>%
  rename(percent_covered = percent_covered, covered_bases = covered_bases, total_bases = total_bases)

# Create output folder if needed
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
  cat("Created output directory:", output_folder, "\n")
}

# Create plots
violin_plot <- plot_breadth_violin(df)
bar_plot <- plot_breadth_bar(df)

violin_file <- file.path(output_folder, paste0(species, "_violin_per_scaffold_percent_covered.png"))
bar_file <- file.path(output_folder, paste0(species, "_barplot_aggregated_percent_covered.png"))

ggsave(violin_file, plot = violin_plot, width = 12, height = 6)
ggsave(bar_file, plot = bar_plot, width = 12, height = 6)

cat("Saved violin plot to:", violin_file, "\n")
cat("Saved bar plot to:", bar_file, "\n")
