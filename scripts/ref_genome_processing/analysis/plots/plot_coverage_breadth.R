library(ggplot2)
library(dplyr)
library(readr)
library(tools)

plot_coverage_breadth_violoin <- function(df_breadth, species) {
  df_breadth$species <- species
  
  plot_breadth <- ggplot(df_breadth, aes(x = factor(species), y = percent_covered)) +
    geom_violin(scale = "width") +
    theme_bw() +
    ylab("Percent Covered") +
    xlab("Species") +
    ggtitle("Distribution of Percent Covered") +
    theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
          legend.text = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          legend.position = "none") +
    scale_fill_viridis_d()

    return(plot_breadth)
}

plot_coverage_breadth <- function(species, filepath, target_folder) {

  # Define output file constants
  PLOT_BREADTH_VIOLIN <- "plot_breadth_violoin.png"
  PLOT_BREADTH_BINS <- "plot_breadth_bins.png"

  # Read the TSV file into a data frame
  df <- read.csv(filepath, header = TRUE)

  # Check if data frame is empty
  if (nrow(df) == 0) {
    stop("WARNING: The input file is empty.")
  }

  # Check if the target folder exists; if not, create it
  if (!dir.exists(target_folder)) {
    dir.create(target_folder, recursive = TRUE)
    cat("Created directory:", target_folder, "\n") 
  }

  # Check and generate plot for coverage breadth violoin
  if (!file.exists(file.path(target_folder, PLOT_BREADTH_VIOLIN))) {  
    plot_coverage_breadth_violoin_var <- plot_coverage_breadth_violoin(df, species)
    save_plot(plot_coverage_breadth_violoin_var, target_folder, PLOT_BREADTH_VIOLIN)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_BREADTH_VIOLIN, "\n")
  }

  # Check and generate plot for coverage breadth bins
  if (!file.exists(file.path(target_folder, PLOT_BREADTH_BINS))) {  
    plot_coverage_breadth_bins_var <- plot_coverage_breadth_bins(df)
    save_plot(plot_coverage_breadth_bins_var, target_folder, PLOT_BREADTH_BINS)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_BREADTH_BINS, "\n")
  }
  
}

save_plot <- function(plot, target_folder, file_name){
  ggsave(file.path(target_folder, file_name), plot = plot, width = 10, height = 6)
  cat("✅ Saved:", file_name, "\n")
}

plot_coverage_breadth_bins <- function(df_breadth) {
  # Bin the scaffolds based on their length (total_bases)
  bins <- c(0, 100000, 250000, 500000, 1000000, 2500000, 5000000, 10000000, 20000000, Inf)  # Updated bins for scaffold length
  bin_labels <- c('0-100k', '100k-250k', '250k-500k', '500k-1M', '1M-2.5M', '2.5M-5M', '5M-10M', '10M-20M', '20M+')

  # Create a new column for the length bin
  df_breadth$length_bin <- cut(df_breadth$total_bases, breaks = bins, labels = bin_labels, right = FALSE)

  # Calculate the average percent_covered, count of scaffolds, and standard deviation for each bin
  avg_coverage_by_bin <- df_breadth %>%
    group_by(length_bin) %>%
    summarise(
      avg_coverage = mean(percent_covered, na.rm = TRUE),
      scaffold_count = n(),
      std_dev = sd(percent_covered, na.rm = TRUE)
    )

  # Create the plot with error bars for standard deviation
  plot <- ggplot(avg_coverage_by_bin, aes(x = length_bin, y = avg_coverage, group = 1)) +
    geom_bar(stat = "identity", fill = "skyblue", color = "black") +
    geom_errorbar(aes(ymin = avg_coverage - std_dev, ymax = avg_coverage + std_dev), 
                  width = 0.2, color = "black") +  # Error bars for standard deviation
    labs(x = "Scaffold Length Bin", y = "Average Percent Covered", 
          title = paste("Average Coverage by Scaffold Length Bin")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

  return(plot)
}

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are passed
if (length(args) < 3) {
  stop("❌ Not enough arguments. Required: species, depth_file, target_folder.")
}

species <- args[1]  # Species (not used in the plot but passed as an argument)
filepath <- args[2]  # Path to the input TSV file
target_folder <- args[3]  # Target folder for saving the plot

cat("🎯 Plotting coverage breadth for species ", species)
cat("Input file: ", filepath)
cat("Target folder: ", target_folder)

plot_coverage_breadth(species, filepath, target_folder)
