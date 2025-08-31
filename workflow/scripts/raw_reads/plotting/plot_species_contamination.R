## Kraken2 contamination boxplot ##

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

#------------------------------------------
# Function: Plot species counts as boxplots
#------------------------------------------
plot_kraken2_barplots <- function(df, species_label, target_folder) {
  base_name <- paste0(species_label, "_kraken2_contamination")
  grouped_path <- file.path(target_folder, paste0(base_name, "_grouped_by_species.png"))
  stacked_path <- file.path(target_folder, paste0(base_name, "_stacked_by_individuum.png"))

  id_col <- "individuum"
  if (!id_col %in% colnames(df)) {
    stop("Column 'individuum' not found in input file.")
  }

  # Convert wide to long format
  df_long <- df %>%
    pivot_longer(
      cols = -all_of(id_col),
      names_to = "species_id",
      values_to = "count"
    ) %>%
    mutate(
      species_id = as.factor(species_id),
      count = replace_na(count, 0)
    )

  #------------------------------------------
  # Grouped Bar Plot: Species on x-axis
  #------------------------------------------
  if (!file.exists(grouped_path)) {
    species_order <- df_long %>%
      group_by(species_id) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      arrange(desc(total)) %>%
      pull(species_id)

    df_long$species_id <- factor(df_long$species_id, levels = species_order)

    p_grouped <- ggplot(df_long, aes(x = species_id, y = count, fill = individuum)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_y_continuous(labels = comma) +
      labs(
        title = paste("Grouped Kraken2 Read Counts by Species -", species_label),
        x = "Species ID",
        y = "Read Count",
        fill = "Individuum"
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    ggsave(grouped_path, plot = p_grouped, width = 12, height = 6, dpi = 300)
    message("Saved grouped barplot to: ", grouped_path)
  } else {
    message("Skipping grouped barplot (already exists): ", grouped_path)
  }

  #------------------------------------------
  # Stacked Barplot: Individual on x-axis
  #------------------------------------------
  if (!file.exists(stacked_path)) {
    individuum_order <- df_long %>%
      group_by(individuum) %>%
      summarise(total = sum(count), .groups = "drop") %>%
      arrange(desc(total)) %>%
      pull(individuum)

    df_long$individuum <- factor(df_long$individuum, levels = individuum_order)

    p_stacked <- ggplot(df_long, aes(x = individuum, y = count, fill = species_id)) +
      geom_bar(stat = "identity") +
      scale_y_continuous(labels = comma) +
      labs(
        title = paste("Stacked Kraken2 Read Counts per Individuum -", species_label),
        x = "Individuum",
        y = "Read Count",
        fill = "Species ID"
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    ggsave(stacked_path, plot = p_stacked, width = 12, height = 6, dpi = 300)
    message("Saved stacked barplot to: ", stacked_path)
  } else {
    message("Skipping stacked barplot (already exists): ", stacked_path)
  }
}


#------------------------------------------
# Master Function: Load data and call plot
#------------------------------------------
plot_kraken2_contamination_summary <- function(species, source_file, target_folder) {
  message("Running plot_kraken2_contamination_summary()")

  if (!file.exists(source_file)) {
    stop("Source file does not exist: ", source_file)
  }

  # Read input file (comma-separated)
  df <- read.csv(source_file, header = TRUE)

  # Check that individuum column exists
  if (!"individuum" %in% colnames(df)) {
    stop("Column 'individuum' not found in input file.")
  }

  # Generate plot
  plot_kraken2_barplots(df, species, target_folder)
}

#------------------------------------------
# Main CLI execution
#------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_kraken2_species_boxplot.R <species> <input_file.csv> <output_folder>")
}

species <- args[1]
input_file <- args[2]
output_folder <- args[3]

plot_kraken2_contamination_summary(species, input_file, output_folder)
