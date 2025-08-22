#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tools)

# ----------------- Plot functions -----------------

# Plot: Minor Allele Frequency (MAF) density
plot_maf_density <- function(df) {
  ggplot(df, aes(x = knownEM)) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    xlab("Minor Allele Frequency (MAF)") +
    ylab("Density") +
    ggtitle("MAF Density Distribution") +
    theme_bw(base_size = 14)
}

# Plot: Histogram of MAFs
plot_maf_histogram <- function(df) {
  ggplot(df, aes(x = knownEM)) +
    geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7, color = "black") +
    xlab("Minor Allele Frequency (MAF)") +
    ylab("Count") +
    ggtitle("Histogram of MAF Values") +
    theme_bw(base_size = 14)
}

# Plot: SNP density per scaffold
plot_snp_density <- function(df) {
  df %>%
    count(chromo) %>%
    ggplot(aes(x = reorder(chromo, -n), y = n, fill = chromo)) +
    geom_col(show.legend = FALSE) +
    xlab("Scaffold") +
    ylab("Number of SNPs") +
    ggtitle("SNP Density per Scaffold") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Plot: SNP density along each scaffold sequence
plot_snp_density_sequence <- function(df) {
  ggplot(df, aes(x = position)) +
    geom_density(fill = "darkorange", alpha = 0.5) +
    facet_wrap(~ chromo, scales = "free_x") +
    xlab("Position on Scaffold") +
    ylab("SNP Density") +
    ggtitle("SNP Density Along Sequences") +
    theme_bw(base_size = 14)
}

# Plot: Reference vs Major allele
plot_ref_vs_major <- function(df) {
  df %>%
    count(ref, major) %>%
    ggplot(aes(x = ref, y = major, fill = n)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c() +
    xlab("Reference Base") +
    ylab("Major Allele") +
    ggtitle("Reference vs. Major Allele") +
    theme_bw(base_size = 14)
}

# ----------------- Utility functions -----------------

# Generate plot only if output file does not already exist
generate_plot_if_needed <- function(file_path, plot_expr) {
  if (!file.exists(file_path)) {
    cat("🛠 Generating plot:", basename(file_path), "\n")
    return(eval(plot_expr))
  } else {
    cat("⏩ Skipping plot (already exists):", basename(file_path), "\n")
    return(NULL)
  }
}

# Save plot only if it's not NULL
save_plot_if_needed <- function(plot, filename, width, height) {
  if (!is.null(plot)) {
    ggsave(filename, plot, width = width, height = height)
    cat("✅ Saved:", filename, "\n")
  }
}

# ----------------- Main -----------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("❌ Usage: Rscript plot_maf_data.R <species> <input.mafs.gz> <output_folder>")
}

species <- args[1]
input_file <- args[2]
output_folder <- args[3]

cat("🎯 Plotting MAF data for ", species, "\n")
cat("Input file: ", input_file, "\n")
cat("Output folder: ", output_folder, "\n")

# Read .mafs.gz file
df <- read_tsv(input_file, comment = "#", col_types = cols())

if (nrow(df) == 0) {
  stop("❌ Input MAF file is empty.")
}

# Explicitly convert 'knownEM' to numeric
df <- df %>%
  mutate(knownEM = as.numeric(knownEM))

  # Optional: Check parsing problems
if (any(is.na(df$knownEM))) {
  cat("⚠️ Some 'knownEM' values could not be parsed. Use problems(df) to debug.\n")
  problems(df)
}

# Create output folder if needed
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
  cat("📁 Created output directory:", output_folder, "\n")
}

# Set file prefix
prefix <- file_path_sans_ext(file_path_sans_ext(basename(input_file)))

# Define output file paths
density_file <- file.path(output_folder, paste0(prefix, "_maf_density.png"))
histogram_file <- file.path(output_folder, paste0(prefix, "_maf_histogram.png"))
snp_density_file <- file.path(output_folder, paste0(prefix, "_snp_density_per_scaffold.png"))
snp_density_seq_file <- file.path(output_folder, paste0(prefix, "_snp_density_along_sequence.png"))
ref_major_file <- file.path(output_folder, paste0(prefix, "_ref_vs_major.png"))

# Generate and save plots (conditionally)
save_plot_if_needed(generate_plot_if_needed(density_file, quote(plot_maf_density(df))), density_file, 8, 5)
save_plot_if_needed(generate_plot_if_needed(histogram_file, quote(plot_maf_histogram(df))), histogram_file, 8, 5)
save_plot_if_needed(generate_plot_if_needed(snp_density_file, quote(plot_snp_density(df))), snp_density_file, 10, 6)
save_plot_if_needed(generate_plot_if_needed(snp_density_seq_file, quote(plot_snp_density_sequence(df))), snp_density_seq_file, 10, 5)
save_plot_if_needed(generate_plot_if_needed(ref_major_file, quote(plot_ref_vs_major(df))), ref_major_file, 6, 5)

cat("🎯 All requested plots are now in:", output_folder, "\n")
