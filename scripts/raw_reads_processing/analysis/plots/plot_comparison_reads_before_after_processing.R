## aDNA sequence length distribution plot ##

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)  # for number formatting

#------------------------------------------
# Function: Plot read count by protocol
#------------------------------------------
plot_by_protocol <- function(df, species, target_folder) {
  file_name <- paste0(species, "_read_count_comparison_protocol.png")
  file_path <- file.path(target_folder, file_name)

  # Skip plotting if file already exists
  if (file.exists(file_path)) {
    message("Skipping protocol plot: already exists at ", file_path)
    return()
  }

  # Summarize counts by protocol
  df_protocol <- df %>%
    mutate(protocol_collapsed = gsub("[12]$", "", protocol)) %>%
    group_by(protocol_collapsed) %>%
    summarise(
      raw_count_absolute = sum(raw_count),
      adapter_removed_count_absolute = sum(adapter_removed_count),
      duplicates_removed_count_absolute = sum(duplicates_removed_count),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(raw_count_absolute, adapter_removed_count_absolute, duplicates_removed_count_absolute),
      names_to = "read_type", values_to = "count"
    ) %>%
    mutate(read_type = factor(
      read_type, levels = c("raw_count_absolute", "adapter_removed_count_absolute", "duplicates_removed_count_absolute")
    ))

  # Generate the bar plot
  p <- ggplot(df_protocol, aes(x = protocol_collapsed, y = count, fill = read_type)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = c(
      "raw_count_absolute" = "#1f77b4",  # Blue
      "adapter_removed_count_absolute" = "#ff7f0e",  # Orange
      "duplicates_removed_count_absolute" = "#2ca02c"  # Green
    )) +
    scale_y_continuous(labels = comma) +
    labs(x = "Protocol", y = "Read Count", fill = "Read Type") +
    theme_bw() +
    theme(panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank())

  # Save plot to file
  ggsave(file_path, plot = p, width = 8, height = 5, dpi = 300)
  message("Saved protocol plot to: ", file_path)
}

#------------------------------------------
# Function: Plot read count by individual
#------------------------------------------
plot_by_individual <- function(df, species, target_folder) {
  file_name <- paste0(species, "_read_count_comparison_individual.png")
  file_path <- file.path(target_folder, file_name)

  # Skip plotting if file already exists
  if (file.exists(file_path)) {
    message("Skipping individual plot: already exists at ", file_path)
    return()
  }

  # Summarize counts by individual
  df_individual <- df %>%
    group_by(individual) %>%
    summarise(
      raw_count_absolute = sum(raw_count),
      adapter_removed_count_absolute = sum(adapter_removed_count),
      duplicates_removed_count_absolute = sum(duplicates_removed_count),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(raw_count_absolute, adapter_removed_count_absolute, duplicates_removed_count_absolute),
      names_to = "read_type", values_to = "count"
    ) %>%
    mutate(read_type = factor(
      read_type, levels = c("raw_count_absolute", "adapter_removed_count_absolute", "duplicates_removed_count_absolute")
    ))

  # Generate the bar plot
  p <- ggplot(df_individual, aes(x = individual, y = count, fill = read_type)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = c(
      "raw_count_absolute" = "#1f77b4",  # Blue
      "adapter_removed_count_absolute" = "#ff7f0e",  # Orange
      "duplicates_removed_count_absolute" = "#2ca02c"  # Green
    )) +
    scale_y_continuous(labels = comma) +
    labs(x = "Individual", y = "Read Count", fill = "Read Type") +
    theme_bw() +
    theme(panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels

  # Save plot to file
  ggsave(file_path, plot = p, width = 8, height = 5, dpi = 300)
  message("Saved individual plot to: ", file_path)
}

#------------------------------------------
# Master Function: Load data and call plot functions
#------------------------------------------
plot_compare_reads_before_after_processing <- function(species, source_file, target_folder) {
  message("Running plot_compare_reads_before_after_processing()")

  # Check if input file exists
  if (!file.exists(source_file)) {
    stop("Source file does not exist: ", source_file)
  }

  # Read the TSV file
  df <- read.table(source_file, header = TRUE)

  # Run individual plot functions
  plot_by_protocol(df, species, target_folder)
  plot_by_individual(df, species, target_folder)
}

#------------------------------------------
# Main CLI execution
#------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Expecting 3 arguments: species, input file, target folder
if (length(args) < 3) {
  stop("Usage: Rscript script_name.R <species> <input_file> <target_folder>")
}

species <- args[1]
input_file <- args[2]
target_folder <- args[3]

plot_compare_reads_before_after_processing(species, input_file, target_folder)
