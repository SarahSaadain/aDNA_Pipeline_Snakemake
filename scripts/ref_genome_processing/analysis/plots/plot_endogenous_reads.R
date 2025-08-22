#load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)  # for number formatting

plot_endogenous_reads <- function(species, source_file, target_folder) {
  
  cat("Executing plot_endogenous_reads", "\n")
  
  # Read the source file
  df <- read.table(
    source_file, 
    sep =",", 
    header = TRUE) # Changed to TRUE to read the header
  
  # Add the Non-Endogenous reads and their percentages to the original dataframe
  df <- df %>%
    mutate(
      Non_Endogenous = TotalReads - MappedReads,  # Changed to use new column names
      Proportion = MappedReads / TotalReads,
      percent_endogenous = Proportion * 100, 
      percent_non_endogenous = (Non_Endogenous / TotalReads) * 100  # Calculate percentage
    )
  
  # Extract the sample name by removing everything from the first dot (.) onward in the 'Filename' column
  df$Individual <- sub("\\..*", "", df$Filename)
  
  # Loop over each row of the dataframe
  for(i in 1:nrow(df)) {
    # Subset the row
    row_data <- df[i, ]

    # Create file name and path for each chart
    file_name <- paste0(row_data$Individual, "_endogenous_reads_pie_chart.png") # Changed to use Filename
    file_path <- file.path(target_folder, file_name)

    if (!file.exists(file_path)) {  # Check if the file already exists
      cat("Generating plot for:", row_data$Individual, "\n")
    } else {
      cat("⏩ Skipping plot (already exists):", row_data$Individual, "\n")
      next  # Skip to the next iteration if the file exists
    }
    
    # Create a long format for the row (pie chart data)
    row_long <- data.frame(
      read_type = c("Endogenous", "Non-Endogenous"),
      count = c(row_data$MappedReads, row_data$Non_Endogenous), # Changed to use new column names
      percent = c(row_data$percent_endogenous, row_data$percent_non_endogenous)
    )
    
    # Create the pie chart
    endogenous_plot_pie_file <- ggplot(row_long, aes(x = "", y = count, fill = read_type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +  # Create the pie chart by converting to polar coordinates
      geom_text(aes(label = paste0(round(percent, 1), "%")), 
                position = position_stack(vjust = 0.5), size = 5, color = "white") +
      scale_fill_manual(values = c("Endogenous" = "#209557",  # Green
                                  "Non-Endogenous" = "#1f5bb4")) + # Blue
      labs(x = NULL, y = NULL, fill = "Read Type", title = paste("Endogenous vs Non-Endogenous Reads:", row_data$Filename)) +  # Include filename in title
      theme_bw() +  # Apply the black-and-white theme
      theme(
        panel.grid = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank(),  # Remove border around the panel
        plot.border = element_blank(),   # Remove border around the entire plot
        axis.text.y = element_blank(),   # Remove y-axis text
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.ticks.y = element_blank(),  # Remove y-axis ticks
        axis.text.x = element_blank()    # Remove x-axis labels if needed
      )
    
    # Save the plot
    ggsave(file_path, plot = endogenous_plot_pie_file, width = 6, height = 6, dpi = 300)

    cat("✅ Saved:", file_name, "\n")
  }

  # Create file name and path for each chart
  file_name <- paste0(species, "_endogenous_reads_bar_chart.png") # Changed to use Filename
  file_path <- file.path(target_folder, file_name)

  if (!file.exists(file_path)) {  # Check if the file already exists
    cat("Generating endogenous reads bar chart for:", species, "\n")
  } else {
    cat("⏩ Skipping plot (already exists): ", file_name, "\n")
    return()  # Skip to the next iteration if the file exists
  }

  # Create the plot
  endogenous_plot_bar <- ggplot(df, aes(x = Individual, y = percent_endogenous, fill = Individual)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    labs(x = "File", y = "Percentage of Endogenous Reads", title = paste("Endogenous Reads")) +
    theme(axis.text.x = element_text(size = 14, angle = 45, vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black")) +
    scale_fill_viridis_d()


  ggsave(file_path, plot = endogenous_plot_bar, width = 6, height = 6, dpi = 300)

  cat("✅ Saved:", file_name, "\n")
}

# FOR TESTING
#plot_endogenous_reads(
# "Bger",
# "/Users/ssaadain/Documents/aDNA/Bger/results/endogenous_reads/Bger_endogenous_reads.csv",
# "/Users/ssaadain/Documents/aDNA/Bger/results/plots/endogenous_reads"
#)


# Ensure you capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are passed
if (length(args) < 3) {
  stop("❌ Not enough arguments. Required: species, depth_file, target_folder.")
}

# Assign the arguments to variables
species <- args[1]
reads_analysis_file <- args[2]
target_folder <- args[3]

cat("🎯 Plotting endogenous reads for species", species, "\n")
cat("Reads analysis file:", reads_analysis_file, "\n")
cat("Target folder:", target_folder, "\n")

plot_endogenous_reads(
  species,
  reads_analysis_file,
  target_folder
)
