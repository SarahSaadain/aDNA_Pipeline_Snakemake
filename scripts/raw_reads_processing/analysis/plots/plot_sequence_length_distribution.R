# Load required libraries
library(ggplot2)
library(tidyr)
library(readr)
library(tools)  # For file path manipulation
library(dplyr)  # For data manipulation
library(scales)  # For number formatting

# Function to generate and save plots
plot_read_length_distribution <- function(species, source_file, output_folder) {

  print("Executing plot_read_length_distribution")
  message("Generating read length distribution plots for species: ", species)

   # Read the TSV data
  df <- read_tsv(source_file, col_types = cols(), show_col_types = FALSE)  # Suppress messages
  
  # Check if required columns exist
  required_cols <- c("reads_file", "read_length", 
                     "read_count_adapter_removed", "read_count_quality_filtered", "read_count_duplicates_removed")
  if (!all(required_cols %in% colnames(df))) {
    stop("Error: The input file must contain the columns: reads_file, read_length, read_count_adapter_removed, read_count_quality_filtered, read_count_duplicates_removed")
  }
  
  # Ensure the output folder exists
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Reshape data from wide to long format
  df_long <- pivot_longer(df, cols = c("read_count_adapter_removed", "read_count_quality_filtered", "read_count_duplicates_removed"),
                          names_to = "Processing_Step", values_to = "Count")

  # Generate a plot for each unique file
  unique_files <- unique(df$reads_file)
  for (file_name in unique_files) {

    message("Generating plot for file: ", file_name)

    # Generate output filename
    output_file <- file.path(output_folder, paste0(file_name, ".png"))

    # Check if the file already exists
    if (file.exists(output_file)) {
      message("File already exists, skipping plot generation for: ", file_name)
      next  # Skip to the next iteration if the file exists
    }

    df_subset <- df_long %>% filter(reads_file == file_name)  # Filter data for this file
    
    p <- ggplot(df_subset, aes(x = read_length, y = Count, color = Processing_Step, group = Processing_Step)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Read Length Distribution:", file_name),
          x = "Read Length",
          y = "Count",
          color = "Processing Step") +
      scale_x_continuous(
        breaks = seq(min(df_subset$read_length, na.rm = TRUE), max(df_subset$read_length, na.rm = TRUE), by = 10),
        labels = label_number()  # This will prevent scientific notation
      ) +
      scale_y_continuous(labels = scales::comma) +  # Format y-axis labels with commas
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels to avoid overlap
      )
    
    # Save the plot
    ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
    
    # Print message
    message("Plot saved to: ", output_file)
  }

  # Generate output filename for the combined plot (same as input file, replacing .tsv with .png)
  combined_output_file <- file.path(output_folder, paste0(file_path_sans_ext(basename(source_file)), ".png"))

  # Check if the combined plot file already exists
  if (file.exists(combined_output_file)) {
    message("Combined plot file already exists, skipping plot generation")
    return()  # Skip to the next iteration if the file exists
  }

  # Aggregate (sum up) counts for the combined plot across all files
  df_combined <- df_long %>%
    group_by(read_length, Processing_Step) %>%
    summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

  # Create a combined plot
  p_combined <- ggplot(df_combined, aes(x = read_length, y = Count, color = Processing_Step, group = Processing_Step)) +
    geom_line(alpha = 0.7) +
    geom_point(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Combined Read Length Distribution Across All Files",
         x = "Read Length",
         y = "Total Count",
         color = "Processing Step") +
    scale_x_continuous(
        breaks = seq(min(df_subset$read_length, na.rm = TRUE), max(df_subset$read_length, na.rm = TRUE), by = 10),
        labels = label_number()  # This will prevent scientific notation
      ) +
      scale_y_continuous(labels = scales::comma) +  # Format y-axis labels with commas+
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels to avoid overlap
      )
  
  # Save the combined plot
  ggsave(combined_output_file, plot = p_combined, width = 10, height = 8, dpi = 300)

  # Print message
  message("Combined plot saved to: ", combined_output_file)
}

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if two arguments are provided
if (length(args) != 3) {
  stop("Not enough arguments. Required: species, depth_file, target_folder.")
}

# Assign the arguments to variables
species <- args[1]
source_file <- args[2]
output_folder <- args[3]

# Run the function
plot_read_length_distribution(species, source_file, output_folder)
