library(ggplot2)
library(stringr)
library(ggpubr)
library(dplyr)
library(purrr)
library(tools)

plot_depth_coverage_violon <- function(df_depth, species) {

  df_depth$species <- species

  # Plot depth
  plot_depth_violin <- ggplot(df_depth, aes(x = factor(species), y = avg_depth)) +
    scale_y_continuous(
      trans = "log10",
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::comma
    ) +
    geom_violin(scale = "width") +
    theme_bw() +
    ylab("Avg. Depth") +
    xlab("Species") +
    ggtitle("Distribution of Average Depth") +
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

    return(plot_depth_violin)
}

# Define the function
plot_depth_coverage <- function(depth, lengthScaffoldRange) {
  
  # Filter the data based on the scaffold length range
  depth_filtered <- depth %>%
    filter(total_bases >= lengthScaffoldRange[1] & total_bases <= lengthScaffoldRange[2])
  
  # Summarize the filtered data by average depth
  depth_filtered_mean <- depth_filtered %>%
    group_by(rounded_avg_depth) %>%
    summarise(nr_scaffolds = n(), .groups = "drop")
  
  # Create the plot
  plot_depth_coverage <- ggplot(depth_filtered_mean, aes(x = rounded_avg_depth, y = nr_scaffolds)) + 
    geom_line(color = "blue", size = 0.3) +  # Line for Mean Depth
    geom_point(color = "blue", size = 0.1) +  # Points for Mean Depth
    labs(title = paste("Depth Coverage of Scaffolds",  
                       " (Length range: ", 
                       format(lengthScaffoldRange[1], big.mark = ",", scientific = FALSE), 
                       " - ", format(lengthScaffoldRange[2], big.mark = ",", scientific = FALSE), ")", sep = ""),
         subtitle = paste(format(nrow(depth_filtered), big.mark = ",", scientific = FALSE), 
                          " of ", format(nrow(depth), big.mark = ",", scientific = FALSE), " Scaffolds", sep = ""),
         x = "Depth", y = "Number of Scaffolds") +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() +
    theme(legend.title = element_blank())

    return(plot_depth_coverage)
}

save_plot <- function(plot, target_folder, file_name){
  ggsave(file.path(target_folder, file_name), plot = plot, width = 10, height = 6)
  cat("✅ Saved:", file_name, "\n")
}

# Define the function
plot_max_depth_coverage <- function(depth, lengthScaffoldRange) {
  
  # Filter based on bin size range
  depth_filtered <- depth[depth$total_bases >= lengthScaffoldRange[1] & depth$total_bases <= lengthScaffoldRange[2], ]
  
  depth_filtered_summary <- depth_filtered %>%
    group_by(rounded_max_depth) %>%
    summarise(nr_scaffolds=n())
  
  # Generate the plot
  depth_plot <- ggplot(depth_filtered_summary, aes(x = rounded_max_depth, y = nr_scaffolds)) +
    geom_line(size = 0.3) +  # Line for trend
    geom_point(size = 0.1) +   # Dots for individual points
    labs(title = paste(
      "Maximal Depth Coverage of Scaffolds",  
      " (Scaffold length from ",  
      format(lengthScaffoldRange[1], big.mark = ",", scientific = FALSE), 
      " - ", 
      format(lengthScaffoldRange[2], big.mark = ",", scientific = FALSE), 
      ")", 
      sep = ""), 
      subtitle = paste( 
        format(nrow(depth_filtered), big.mark = ",", scientific = FALSE), 
        " of ",
        format(nrow(depth), big.mark = ",", scientific = FALSE), 
        " Scaffolds",
        sep = ""), 
      x = "Maximum depth", y = "Number of Scaffolds") +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw() +
    theme(legend.position = "none")
  
  # Return the plot
  return(depth_plot)
}

plot_coverage_depth <- function(species, depth_breath_analysis_file_path, target_folder){

  # Define constants for file names
  PLOT_DEPTH_COVERAGE_ALL <- "plot_DepthCoverageOfAllScaffolds.png"
  PLOT_DEPTH_COVERAGE_SUPER_SMALL <- "plot_DepthCoverageOfSuperSmallScaffolds.png"
  PLOT_DEPTH_COVERAGE_SMALL <- "plot_DepthCoverageOfSmallScaffolds.png"
  PLOT_DEPTH_COVERAGE_MEDIUM <- "plot_DepthCoverageOfMediumScaffolds.png"
  PLOT_DEPTH_COVERAGE_LARGE <- "plot_DepthCoverageOfLargeScaffolds.png"

  PLOT_MAX_DEPTH_COVERAGE_ALL <- "plot_MaxDepthCoverageOfAllScaffolds.png"
  PLOT_MAX_DEPTH_COVERAGE_SUPER_SMALL <- "plot_MaxDepthCoverageOfSuperSmallScaffolds.png"
  PLOT_MAX_DEPTH_COVERAGE_SMALL <- "plot_MaxDepthCoverageOfSmallScaffolds.png"
  PLOT_MAX_DEPTH_COVERAGE_MEDIUM <- "plot_MaxDepthCoverageOfMediumScaffolds.png"
  PLOT_MAX_DEPTH_COVERAGE_LARGE <- "plot_MaxDepthCoverageOfLargeScaffolds.png"

  PLOT_DEPTH_COVERAGE_VIOLIN = "plot_depthCoverage_violin.png"
  
  cat("Executing plot_coverage_depth", "\n")
  
  # depth coverage for each scaffold
  df_depth <- read.csv(depth_breath_analysis_file_path, header = TRUE)

  if (nrow(df_depth) == 0) {
    stop("WARNING: Input file is empty.")
  }

    # Check if the target folder exists; if not, create it
  if (!dir.exists(target_folder)) {
    dir.create(target_folder, recursive = TRUE)
    cat("Created directory:", target_folder, "\n")
  }
  
  # Calculate the mean depth for each scaffold
  df_depth <- df_depth %>%
    mutate(
      rounded_avg_depth = round(avg_depth),
      rounded_max_depth = round(max_depth)
    )
  
  max_length_of_scaffold <- max(df_depth$total_bases, na.rm = TRUE)

  # Check and generate/save plot for depth coverage
  if (!file.exists(file.path(target_folder, PLOT_DEPTH_COVERAGE_ALL))) {
    plot_DepthCoverageOfAllScaffolds <- plot_depth_coverage(df_depth, c(0, max_length_of_scaffold))
    save_plot(plot_DepthCoverageOfAllScaffolds, target_folder, PLOT_DEPTH_COVERAGE_ALL)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_DEPTH_COVERAGE_ALL, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_DEPTH_COVERAGE_SUPER_SMALL))) {
    plot_DepthCoverageOfSuperSmallScaffolds <- plot_depth_coverage(df_depth, c(0, 1000))
    save_plot(plot_DepthCoverageOfSuperSmallScaffolds, target_folder, PLOT_DEPTH_COVERAGE_SUPER_SMALL)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_DEPTH_COVERAGE_SUPER_SMALL, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_DEPTH_COVERAGE_SMALL))) {
    plot_DepthCoverageOfSmallScaffolds <- plot_depth_coverage(df_depth, c(1001, 10000))
    save_plot(plot_DepthCoverageOfSmallScaffolds, target_folder, PLOT_DEPTH_COVERAGE_SMALL)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_DEPTH_COVERAGE_SMALL, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_DEPTH_COVERAGE_MEDIUM))) {
    plot_DepthCoverageOfMediumScaffolds <- plot_depth_coverage(df_depth, c(10001, 100000))
    save_plot(plot_DepthCoverageOfMediumScaffolds, target_folder, PLOT_DEPTH_COVERAGE_MEDIUM)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_DEPTH_COVERAGE_MEDIUM, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_DEPTH_COVERAGE_LARGE))) {
    plot_DepthCoverageOfLargeScaffolds <- plot_depth_coverage(df_depth, c(100001, max_length_of_scaffold))
    save_plot(plot_DepthCoverageOfLargeScaffolds, target_folder, PLOT_DEPTH_COVERAGE_LARGE)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_DEPTH_COVERAGE_LARGE, "\n")
  }

  # Check and generate/save plot for max depth coverage
  if (!file.exists(file.path(target_folder, PLOT_MAX_DEPTH_COVERAGE_ALL))) {
    plot_MaxDepthCoverageOfAllScaffolds <- plot_max_depth_coverage(df_depth, c(0, max_length_of_scaffold))
    save_plot(plot_MaxDepthCoverageOfAllScaffolds, target_folder, PLOT_MAX_DEPTH_COVERAGE_ALL)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_MAX_DEPTH_COVERAGE_ALL, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_MAX_DEPTH_COVERAGE_SUPER_SMALL))) {
    plot_MaxDepthCoverageOfSuperSmallScaffolds <- plot_max_depth_coverage(df_depth, c(0, 1000))
    save_plot(plot_MaxDepthCoverageOfSuperSmallScaffolds, target_folder, PLOT_MAX_DEPTH_COVERAGE_SUPER_SMALL)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_MAX_DEPTH_COVERAGE_SUPER_SMALL, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_MAX_DEPTH_COVERAGE_SMALL))) {
    plot_MaxDepthCoverageOfSmallScaffolds <- plot_max_depth_coverage(df_depth, c(1001, 10000))
    save_plot(plot_MaxDepthCoverageOfSmallScaffolds, target_folder, PLOT_MAX_DEPTH_COVERAGE_SMALL)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_MAX_DEPTH_COVERAGE_SMALL, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_MAX_DEPTH_COVERAGE_MEDIUM))) {
    plot_MaxDepthCoverageOfMediumScaffolds <- plot_max_depth_coverage(df_depth, c(10001, 100000))
    save_plot(plot_MaxDepthCoverageOfMediumScaffolds, target_folder, PLOT_MAX_DEPTH_COVERAGE_MEDIUM)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_MAX_DEPTH_COVERAGE_MEDIUM, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_MAX_DEPTH_COVERAGE_LARGE))) {
    plot_MaxDepthCoverageOfLargeScaffolds <- plot_max_depth_coverage(df_depth, c(100001, max_length_of_scaffold))
    save_plot(plot_MaxDepthCoverageOfLargeScaffolds, target_folder, PLOT_MAX_DEPTH_COVERAGE_LARGE)
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_MAX_DEPTH_COVERAGE_LARGE, "\n")
  }

  if (!file.exists(file.path(target_folder, PLOT_DEPTH_COVERAGE_VIOLIN))) {
    plot_depth_violin <- plot_depth_coverage_violon(df_depth, species)
    save_plot(plot_depth_violin, target_folder, PLOT_DEPTH_COVERAGE_VIOLIN)
    
  } else {
    cat("⏩ Skipping plot (already exists):", PLOT_DEPTH_COVERAGE_VIOLIN, "\n")
  }
  
}

# FOR TESTING
# plot_coverage_depth(
#   "Bger",
#   "/Users/ssaadain/Documents/aDNA/Bger/results/qualitycontrol/depth_breadth/C1.fastq_GCA_000762945.2_Bger_2.0_genomic_analysis.tsv",
#   "/Users/ssaadain/Documents/aDNA/Bger/results/plots/depth/C1.fastq_GCA_000762945.2_Bger_2.0_genomic"
# )

# Ensure you capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are passed
if (length(args) < 3) {
  stop("❌ Not enough arguments. Required: species, depth_file, target_folder.")
}

# Assign the arguments to variables
species <- args[1]
depth_file <- args[2]
target_folder <- args[3]

cat("🎯 Plotting coverage depth for species", species, "\n")
cat("Depth file:", depth_file, "\n")
cat("Target folder:", target_folder, "\n")

plot_coverage_depth(
  species,
  depth_file,
  target_folder
)



