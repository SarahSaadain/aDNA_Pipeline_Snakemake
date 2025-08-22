#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tools)

# ----------------- Plot functions -----------------

# Plot: Minor Allele Frequency (MAF) density comparison across individuals
plot_maf_density_comparison <- function(df) {
  ggplot(df, aes(x = knownEM, color = individual, fill = individual)) +
    geom_density(alpha = 0.3) +
    xlab("Minor Allele Frequency (MAF)") +
    ylab("Density") +
    ggtitle("MAF Density Distribution Across Individuals") +
    theme_bw(base_size = 14)
}

# Plot: Histogram of MAFs comparison across individuals
plot_maf_histogram_comparison <- function(df) {
  ggplot(df, aes(x = knownEM, fill = individual)) +
    geom_histogram(bins = 30, alpha = 0.5, position = "identity", color = "black") +
    xlab("Minor Allele Frequency (MAF)") +
    ylab("Count") +
    ggtitle("Histogram of MAF Values Across Individuals") +
    theme_bw(base_size = 14)
}

# Plot: SNP density per scaffold comparison (faceted by individual)
plot_snp_density_comparison <- function(df) {
  df %>%
    count(individual, chromo) %>%
    ggplot(aes(x = reorder(chromo, -n), y = n, fill = individual)) +
    geom_col(show.legend = TRUE, position = "dodge") +
    xlab("Scaffold") +
    ylab("Number of SNPs") +
    ggtitle("SNP Density per Scaffold Across Individuals") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Plot: SNP density along each scaffold sequence (faceted by chromo and individual)
plot_snp_density_sequence_comparison <- function(df) {
  ggplot(df, aes(x = position, color = individual, fill = individual)) +
    geom_density(alpha = 0.3) +
    facet_grid(chromo ~ individual, scales = "free_x") +
    xlab("Position on Scaffold") +
    ylab("SNP Density") +
    ggtitle("SNP Density Along Sequences Across Individuals") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Plot: Total SNP count per individual (bar plot)
plot_total_snp_count_per_individual <- function(df) {
  df %>%
    count(individual) %>%
    ggplot(aes(x = reorder(individual, -n), y = n, fill = individual)) +
    geom_col(show.legend = FALSE) +
    xlab("Individual") +
    ylab("Total Number of SNPs") +
    ggtitle("Total SNP Count per Individual") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Note: The ref vs major allele heatmap may be less meaningful for multi-individual comparisons
# So we remove or optionally create separate heatmaps per individual if desired.

# ----------------- Utility functions -----------------

generate_plot_if_needed <- function(file_path, plot_expr) {
  if (!file.exists(file_path)) {
    cat("🛠 Generating plot:", basename(file_path), "\n")
    return(eval(plot_expr))
  } else {
    cat("⏩ Skipping plot (already exists):", basename(file_path), "\n")
    return(NULL)
  }
}

save_plot_if_needed <- function(plot, filename, width, height) {
  if (!is.null(plot)) {
    ggsave(filename, plot, width = width, height = height)
    cat("✅ Saved:", filename, "\n")
  }
}

# ----------------- Main -----------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("❌ Usage: Rscript plot_maf_data_multi.R <species> <input_folder> <output_folder>")
}

species <- args[1]
input_folder <- args[2]
output_folder <- args[3]

cat("🎯 Plotting MAF data for", species, "\n")
cat("Input folder:", input_folder, "\n")
cat("Output folder:", output_folder, "\n")

# List all .mafs.gz files in the input folder
maf_files <- list.files(input_folder, pattern = "\\.mafs\\.gz$", full.names = TRUE)

if (length(maf_files) == 0) {
  stop("❌ No .mafs.gz files found in input folder.")
}

# Read all files and add an 'individual' column based on filename
df_list <- lapply(maf_files, function(file) {
  cat("📥 Reading file:", basename(file), "\n")
  df <- read_tsv(file, comment = "#", col_types = cols())
  if (nrow(df) == 0) {
    warning(paste("⚠️ File", basename(file), "is empty. Skipping."))
    return(NULL)
  }
  # take part before the first . of filename  as individual
  df$individual <- strsplit(basename(file), "[.]")[[1]][1]
  return(df)
})

df <- bind_rows(df_list)

if (nrow(df) == 0) {
  stop("❌ No data loaded from files.")
}

# Create output folder if needed
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
  cat("📁 Created output directory:", output_folder, "\n")
}

density_file <- file.path(output_folder, paste0(species, "_maf_density_comparison.png"))
histogram_file <- file.path(output_folder, paste0(species, "_maf_histogram_comparison.png"))
snp_density_file <- file.path(output_folder, paste0(species, "_snp_density_per_scaffold_comparison.png"))
snp_density_seq_file <- file.path(output_folder, paste0(species, "_snp_density_along_sequence_comparison.png"))
total_snp_count_file <- file.path(output_folder, paste0(species, "_total_snp_count_per_individual.png"))

# Generate and save plots (conditionally)
save_plot_if_needed(generate_plot_if_needed(density_file, quote(plot_maf_density_comparison(df))), density_file, 8, 5)
save_plot_if_needed(generate_plot_if_needed(histogram_file, quote(plot_maf_histogram_comparison(df))), histogram_file, 8, 5)
save_plot_if_needed(generate_plot_if_needed(snp_density_file, quote(plot_snp_density_comparison(df))), snp_density_file, 12, 6)
#save_plot_if_needed(generate_plot_if_needed(snp_density_seq_file, quote(plot_snp_density_sequence_comparison(df))), snp_density_seq_file, 12, 8)
save_plot_if_needed(generate_plot_if_needed(total_snp_count_file, quote(plot_total_snp_count_per_individual(df))), total_snp_count_file, 7, 5)

cat("🎯 All requested comparison plots are now in:", output_folder, "\n")
