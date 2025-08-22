library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(yaml) # Load the yaml library

process_and_plot_before_after <- function(analysis_files, output_folder, species_names, comparison_name) {
  list_of_analysis_dataframes <- list()
  
  for (i in 1:length(analysis_files)) {
    filepath <- analysis_files[i]

    print(paste("Processing file:", filepath))

        # Check if the file does not exists exit
    if (!file.exists(filepath)) {
      message(paste("File not found:", filepath))
      return(NULL)  # exit the function early without an error
    } 
    
    df <- read.table(
      filepath,
      sep = "\t",
      header = TRUE,
      colClasses = c(individual = "character", protocol = "character")
    )
    df$species_id <- names(species_names)[i] #get the species ID
    df$species <- species_names[[i]]  # Get species long name from the provided list
    list_of_analysis_dataframes[[df$species_id[1]]] <- df # use species id to store
   
  }
  
  df_before_after <- bind_rows(list_of_analysis_dataframes)

  #check if the species column exists
  if (!("species" %in% colnames(df_before_after))){
    stop("Error: 'species' column not found in the dataframes.  Check the input files.")
  }

  # Create a factor for the 'species' column, defining the order of the levels.
  if (!is.null(species_names)) {
    desired_order = unname(species_names) # use the long names for ordering
    df_before_after$species <- factor(df_before_after$species, levels = desired_order)
  }
  
  df_species <- df_before_after %>%
    group_by(species) %>%
    summarise(raw_count_absolute = sum(raw_count),
              adapter_removed_count_absolute = sum(adapter_removed_count),
              duplicates_removed_count_absolute = sum(duplicates_removed_count),
              .groups = "drop")
  
  df_long_species <- df_species %>%
    pivot_longer(
      cols = c(raw_count_absolute,
               adapter_removed_count_absolute,
               duplicates_removed_count_absolute),
      names_to = "read_type", values_to = "count"
    )
  
  df_long_species <- df_long_species %>%
    mutate(
      read_type = factor(
        read_type, levels = c("raw_count_absolute",
                              "adapter_removed_count_absolute",
                              "duplicates_removed_count_absolute")))

  # Create a factor for the 'species' column, defining the order of the levels.
  if (!is.null(species_names)) {
    desired_order = unname(species_names) # use the long names for ordering
    df_long_species$species <- factor(df_long_species$species, levels = desired_order)
  }                            
  
  plot_before_after <- ggplot(df_long_species, aes(x = species, y = count, fill = read_type)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    scale_fill_manual(
      values = c("raw_count_absolute" = "darkblue",
                 "adapter_removed_count_absolute" = "#007FFF",
                 "duplicates_removed_count_absolute" = "#00FFEF"),
      labels = c("raw_reads",
                 "reads_after_adapterremoval",
                 "reads_after_duplicationremoval")
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Species", y = "Read Count", fill = "Read Type", title = "Read count raw vs after processing") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 18, angle = 45, vjust = 1, hjust = 1),
          legend.text = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"))
  
  ggsave(file.path(output_folder, paste0(comparison_name, "_plot_before_after.png")), plot_before_after, width = 12, height = 8, dpi = 300)
}



# Command line argument handling
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript script_name.R <output_folder> <config_file> <root_folder>")
}

adna_project_folder_path <- args[1]
adna_config_file <- args[2] #  argument for the config file
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

# Loop over each species comparison defined in the config file.
# Each comparison involves a group of species to be analyzed together.
for (comparison_name in names(config$compare_species)) {
  # Extract the data for the current comparison.
  # This 'comparison_data' (referred to as '_data' in the original code)
  comparison_data <- config$compare_species[[comparison_name]]
  
  # Initialize empty lists to store the file paths and species names.
  # Using lists first allows for flexible construction before converting to named vectors.
  analysis_files_list <- list()
  species_names_list <- list()
  
  # Iterate through each entry within the current comparison.
  for (comparison_entry_key in names(comparison_data)) {
    current_entry_details <- comparison_data[[comparison_entry_key]]
    
    # Determine the actual species_id based on the new rule:
    # If 'species_id' is explicitly provided in the current entry and is not empty, use it.
    # Otherwise, default to using the 'comparison_entry_key' itself as the species_id.
    actual_species_id <- if (!is.null(current_entry_details$species_id) && current_entry_details$species_id != "") {
      current_entry_details$species_id
    } else {
      comparison_entry_key
    }
    
    # Retrieve the folder name for the determined 'actual_species_id' from the main config.
    species_folder <- config$species[[actual_species_id]]$folder_name
    
    # Extract the base name (filename without extension) of the reference genome file
    # for the current comparison entry.
    ref_genome_name <- tools::file_path_sans_ext(basename(current_entry_details$reference_genome))
    
    # Construct the full file path to the extended coverage analysis file for this entry.
    # The 'species_id' used in the filename is the 'actual_species_id' determined above.
    file_path <- file.path(
      adna_project_folder_path,
      species_folder,
      "results",
      "qualitycontrol",
      "processed_reads",
      paste0(actual_species_id, "_reads_processing_result.tsv")
    )
    # Store the constructed file path, using the original 'comparison_entry_key'
    # as the name for easy reference later.
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
  # This makes them suitable for functions that expect named vectors.
  analysis_files <- unlist(analysis_files_list)
  species_names <- unlist(species_names_list)
  
  # Check if the list of analysis files is empty after processing all entries.
  # If so, issue a warning and skip to the next comparison.
  if (length(analysis_files) == 0) {
    warning(paste("No analysis files found for comparison:", comparison_name))
    next  # Continue to the next comparison in the loop
  }

  # Log information about the species analyzed to the console for debugging/monitoring.
  cat("Species analyzed for comparison:", comparison_name, "\n")
  cat("Species IDs (as used in config lookup):", names(species_names), "\n") # These are the keys from the comparison data
  cat("Species long names (for display):", species_names, "\n")
  cat("Analysis files:", analysis_files, "\n")
  
  # Log the reference genomes used for each entry in the comparison.
  # This uses the original 'comparison_entry_key' to access the 'reference_genome' from 'comparison_data'.
  cat("Reference genomes:", sapply(names(comparison_data), function(comparison_entry_key) {
    comparison_data[[comparison_entry_key]]$reference_genome
  }), "\n")
  
  process_and_plot_before_after(
    analysis_files,   # Vector of file paths for this comparison
    output_folder_path_for_plots,    # Destination folder for saving the plot
    species_names,    # Named vector of species long names
    comparison_name   # Name of the current comparison, for labeling
  )

  print(paste("Depth and breadth plot saved for comparison:", comparison_name))
}


