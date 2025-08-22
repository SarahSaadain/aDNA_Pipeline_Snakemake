library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(stringr)
library(scales)
library(yaml) # Load the yaml library

process_and_plot_endogenous_reads <- function(analysis_files, output_folder, species_names, comparison_name) {
  # Initialize an empty list to store data frames
  list_of_analysis_dataframes <- list()
  
  # Iterate over each analysis file
  for (i in 1:length(analysis_files)) {
    
    file_path <- analysis_files[i]

    print(paste("Processing file:", file_path))
    
    # Check if the file does not exists exit
    if (!file.exists(file_path)) {
      message(paste("File not found:", file_path))
      return(NULL)  # exit the function early without an error
    } 

    # Read the CSV file, now with header=TRUE
    df <- read.csv(file_path, header = TRUE)
    
    # Assign column names to the expected names
    col_names <- c("sample", "endogenous", "total", "percent_endogenous")

    # Check if the file contains the expected columns.  If not, error.
    if (!all(c("Filename", "MappedReads", "TotalReads", "Proportion") %in% colnames(df))){
      stop(paste("File", file_path, "does not contain the expected columns: Filename, MappedReads, TotalReads, Proportion. Please check the input data."))
    }
    
    # Rename the columns to standard names
    df <- df %>%
      rename(
        sample = Filename,
        endogenous = MappedReads,
        total = TotalReads,
        percent_endogenous = Proportion
      )
    
    df$species_id <- names(species_names)[i]  # Get species ID
    df$species <- species_names[[i]] #get species long name
    
    list_of_analysis_dataframes[[df$species_id[1]]] <- df
  }
  
  # Combine all data frames into one data frame
  df_combined <- bind_rows(list_of_analysis_dataframes)
  
  # Create a new column combining species and sample name for x-axis labels.  Use species long name.
  df_combined$label <- paste(df_combined$species, df_combined$sample, sep = "_")
  
  # Convert percent_endogenous to percentage
  df_combined$percent_endogenous <- df_combined$percent_endogenous * 100
  
  # Create a factor for the 'species' column, defining the order of the levels.
  if (!is.null(species_names)) {
    desired_order = unname(species_names)  # Use long names for ordering
    df_combined$species <- factor(df_combined$species, levels = desired_order)
  }
  
  # Reorder labels based on desired species order
  df_combined$label <- factor(df_combined$label, levels = df_combined$label[order(match(df_combined$species, desired_order))])
  
  # Create the plot
  endogenous_plot <- ggplot(df_combined, aes(x = species, y = percent_endogenous, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw(base_size = 16) +
    labs(x = "Species", y = "Percentage of Endogenous Reads", title = "Endogenous Reads Across Species") +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) +
    scale_fill_viridis_d()
  
  # Print the plot
  #print(endogenous_plot)
  
  # Save the plot to the specified output folder
  ggsave(file.path(output_folder, paste0(comparison_name, "_plot_endogenous_reads.png")), endogenous_plot, width = 12, height = 8, dpi = 300)
}

# Command line argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript script_name.R <root_folder> <config_file> <output_folder>")
}

adna_project_folder_path <- args[1]
adna_config_file <- args[2]
output_folder_path_for_plots <- args[3]

# Create the output directory if it doesn't exist
if (!dir.exists(output_folder_path_for_plots)) {
  dir.create(output_folder_path_for_plots, recursive = TRUE)
}

# Check if the root folder exists
if (!dir.exists(adna_project_folder_path)) {
  stop(paste("Root folder does not exist:", adna_project_folder_path))
}

# Check if the config file exists
if (!file.exists(adna_config_file)) {
  stop(paste("Config file does not exist:", adna_config_file))
}

# Read the config file
config <- yaml.load_file(adna_config_file) # Load the config file

# Check if any relevant configuration is present
if (is.null(config$compare_species) || length(config$compare_species) == 0) {
  stop("No comparisons found in the config file.")
}

# Loop over each species comparison defined in the config's 'compare_species' section.
# Each entry defines a group of species to compare for endogenous reads.
for (comparison_name in names(config$compare_species)) {
  
  # Retrieve the species data for this comparison.
  # Each item includes the species ID and its associated reference genome.
  comparison_data <- config$compare_species[[comparison_name]]
  
  # Initialize empty lists to store the file paths and species names.
  # Using lists first allows for flexible construction before converting to named vectors.
  analysis_files_list <- list()
  species_names_list <- list()
  
  # Iterate through each entry within the current comparison.
  for (comparison_entry_key in names(comparison_data)) {
    current_entry_details <- comparison_data[[comparison_entry_key]]
    
    # Determine the actual species_id to be used for looking up details in config$species.
    # It first checks if 'species_id' is explicitly provided within the current entry's details.
    # If not, it defaults to using the 'comparison_entry_key' itself as the species_id.
    actual_species_id <- if (!is.null(current_entry_details$species_id) && current_entry_details$species_id != "") {
      current_entry_details$species_id
    } else {
      comparison_entry_key
    }
    
    # Get the folder name for the determined 'actual_species_id' from the main config$species section.
    species_folder <- config$species[[actual_species_id]]$folder_name
    
    # Extract the base name (filename without extension) of the reference genome file
    # specified for the current comparison entry.
    ref_genome_name <- tools::file_path_sans_ext(basename(current_entry_details$reference_genome))
    
    # Build the full file path to the endogenous reads CSV file for this species.
    # The 'species_id' used in the filename part of the path is the 'actual_species_id'.
    file_path <- file.path(
      adna_project_folder_path,
      species_folder,
      "results",
      ref_genome_name,
      "endogenous_reads",
      paste0(actual_species_id, "_endogenous_reads.csv")
    )
    # Store the constructed file path in the list, using the original 'comparison_entry_key'
    # as the name. This key is what defines the relationship within the comparison.
    analysis_files_list[[comparison_entry_key]] <- file_path
    
    # Get the user-friendly species name from the main config$species section.
    # This name will be used as a base for generating a unique plot label.
    display_name_from_config <- config$species[[actual_species_id]]$name
    if (is.null(display_name_from_config) || display_name_from_config == "") {
      display_name_from_config <- actual_species_id # Fallback if no display name is found
    }
    
    # Determine the label for plotting (value in species_names_list).
    # This logic aims to create a unique and informative label to prevent "duplicated factor level" errors.
    if (!is.null(current_entry_details$species_id) && current_entry_details$species_id != "") {
        # If species_id was explicitly given, concatenate comparison_entry_key with the display name.
        # This makes the label unique if multiple entries refer to the same actual_species_id.
        species_names_list[[comparison_entry_key]] <- paste0(comparison_entry_key, " ", display_name_from_config)
    } else {
        # If species_id was NOT explicitly given, use the display_name_from_config directly.
        species_names_list[[comparison_entry_key]] <- display_name_from_config
    }
  }
  
  # Convert the lists of file paths and species names into named vectors.
  # This makes them suitable for functions that expect named vectors (e.g., for plotting).
  analysis_files <- unlist(analysis_files_list)
  species_names <- unlist(species_names_list)
  
  # If no analysis files were generated (empty list), issue a warning and skip this comparison.
  if (length(analysis_files) == 0) {
    warning(paste("No analysis files found for comparison:", comparison_name))
    next  # Proceed to the next comparison
  }

  # Information about the species analyzed to console
  cat("Species analyzed for comparison:", comparison_name, "\n")
  # Note: names(species_names) will correspond to the 'comparison_entry_key'
  cat("Comparison Entry Keys:", names(species_names), "\n")
  cat("Species Display Names:", species_names, "\n")
  cat("Analysis files:", analysis_files, "\n")
  # Log the reference genomes used for each entry in the comparison.
  # This uses the original 'comparison_entry_key' to access the 'reference_genome' from 'comparison_data'.
  cat("Reference genomes:", sapply(names(comparison_data), function(comparison_entry_key) {
    comparison_data[[comparison_entry_key]]$reference_genome
  }), "\n")
  
  # Call the function
  process_and_plot_endogenous_reads(
    analysis_files,   # Vector of file paths for this comparison
    output_folder_path_for_plots,    # Destination folder for saving the plot
    species_names,    # Named vector of species long names
    comparison_name   # Name of the current comparison, for labeling
  )

  print(paste("Endogenous reads plot saved for comparison:", comparison_name))
}
