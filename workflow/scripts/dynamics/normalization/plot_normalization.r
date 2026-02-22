#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(optparse)
  library(data.table)
})

# Define command-line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input TSV file [required]", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="output_plots",
              help="Output directory [default=%default]", metavar="DIR"),
  make_option(c("--order"), type="character", default=NULL,
              help="Comma-separated list of individual order (e.g., 'Bger04,Bger05,Bger06')", 
              metavar="STRING"),
  make_option(c("--names"), type="character", default=NULL,
              help="Comma-separated list of display names matching order (optional)", 
              metavar="STRING")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Generate barplots for SCG and TE data")
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input file is required (--input)", call.=FALSE)
}

if (!file.exists(opt$input)) {
  stop(paste("Input file not found:", opt$input), call.=FALSE)
}

# Read the TSV file
# Robust TSV reader (handles CRLF, missing fields, large files)
data_dt <- fread(
  opt$input,
  sep = "\t",
  header = TRUE,
  na.strings = c("NA", ""),
  data.table = FALSE,
  check.names = FALSE
)

# Move first column to rownames
rownames(data_dt) <- data_dt[[1]]
data <- data_dt[, -1, drop = FALSE]

# Determine individual order
if (!is.null(opt$order)) {
  individual_order <- unlist(strsplit(opt$order, ","))
  individual_order <- trimws(individual_order)
} else {
  individual_order <- colnames(data)
}

# Check that all individuals in order exist in data
missing_individuals <- setdiff(individual_order, colnames(data))
if (length(missing_individuals) > 0) {
  stop(paste("Individuals not found in data:", paste(missing_individuals, collapse=", ")), 
       call.=FALSE)
}

# Determine display names (labels) for plotting
if (!is.null(opt$names)) {
  # Use user-provided names as labels exactly
  display_names_vec <- unlist(strsplit(opt$names, ","))
  display_names_vec <- trimws(display_names_vec)
  
  if (length(display_names_vec) != length(individual_order)) {
    stop(paste0("Number of display names (", length(display_names_vec), 
                ") must match number of individuals (", length(individual_order), ")"))
  }
  individual_labels <- setNames(display_names_vec, individual_order)
} else {
  # If no names provided, use the order IDs themselves exactly
  individual_labels <- setNames(individual_order, individual_order)
}

# Separate SCG and TE data
scg_data <- data[grepl("_scg$", rownames(data)), ]
te_data <- data[grepl("_fle$", rownames(data)), ]

cat("Data summary:\n")
cat("  SCG entries:", nrow(scg_data), "\n")
cat("  FLE entries:", nrow(te_data), "\n")
cat("  Individuals:", length(individual_order), "\n\n")

# Function to extract TE family from row names
extract_te_family <- function(rowname) {
  if (grepl("#", rowname)) {
    # Extract the full string, e.g., "LTR/Pao" or just "LTR"
    full_type <- sub(".*#(.*)_fle$", "\\1", rowname)
    
    # If it's just a top-level class (no slash), label it as unclassified
    if (!grepl("/", full_type)) {
      return(paste0(full_type, "/Unclassified"))
    }
    return(full_type)
  }
  return("Unknown")
}

# Function to extract top-level TE class (LINE, DNA, LTR, etc.)
extract_te_class <- function(family) {
  # Extract text before first /
  if (grepl("/", family)) {
    class <- sub("^([^/]+)/.*$", "\\1", family)
    return(class)
  }
  return(family)
}

# Add family and class columns to TE data
te_families <- sapply(rownames(te_data), extract_te_family)
te_data$family <- te_families
te_data$class <- sapply(te_families, extract_te_class)

# Calculate mean values per individual (excluding NAs)
calculate_means <- function(df, exclude_cols = NULL) {
  if (!is.null(exclude_cols)) {
    df_numeric <- df[, !colnames(df) %in% exclude_cols]
  } else {
    df_numeric <- df
  }
  
  means <- colMeans(df_numeric, na.rm = TRUE)
  return(data.frame(
    individual = names(means),
    mean_value = means,
    stringsAsFactors = FALSE
  ))
}

# Calculate means for SCG data
scg_means <- calculate_means(scg_data)
scg_means$type <- "SCG"
scg_means$category <- "SCG"

# Calculate means for all TEs
te_means_all <- calculate_means(te_data, exclude_cols = c("family", "class"))
te_means_all$type <- "All_TEs"
te_means_all$category <- "TE"

# Calculate means per TE family (specific families like DNA/TcMar-Tc1, LINE/L2, etc.)
te_families_unique <- unique(te_data$family)
te_means_by_family <- lapply(te_families_unique, function(fam) {
  family_data <- te_data[te_data$family == fam, ]
  family_data_clean <- family_data[, !colnames(family_data) %in% c("family", "class")]
  means <- calculate_means(family_data_clean)
  means$type <- fam
  means$category <- "TE_family"
  return(means)
})

# Calculate means per TE class (LINE, DNA, LTR, etc.)
te_classes_unique <- unique(te_data$class)
te_means_by_class <- lapply(te_classes_unique, function(cls) {
  class_data <- te_data[te_data$class == cls, ]
  class_data_clean <- class_data[, !colnames(class_data) %in% c("family", "class")]
  means <- calculate_means(class_data_clean)
  means$type <- cls
  means$category <- "TE_class"
  return(means)
})

# Combine all data
all_means <- rbind(
  scg_means, 
  te_means_all, 
  do.call(rbind, te_means_by_family),
  do.call(rbind, te_means_by_class)
)

# Apply custom order and labels - FIX: Make both individual AND display_label factors with correct order
all_means$individual <- factor(all_means$individual, levels = individual_order)
all_means$display_label <- individual_labels[as.character(all_means$individual)]
# Make display_label a factor with the same order to preserve sorting in plots
all_means$display_label <- factor(all_means$display_label, 
                                   levels = individual_labels[individual_order])

# Create output directory
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)

# Function to create barplot
# --- Updated Barplot Function ---

create_barplot <- function(data, title, name_only, folder, prefix = "", keep_all_levels = TRUE) {
  
  # Skip if zero coverage or no data
#   if (nrow(data) == 0 || max(data$mean_value, na.rm = TRUE) == 0) {
#     cat("  Skipping zero-coverage data:", title, "\n")
#     return(NULL)
#   }
  
  # 1. Clean the TE name for the filesystem (e.g., "DNA/Tc1" -> "DNA_Tc1")
  safe_name <- gsub("[/?*:|\"<>\\\\]", "_", name_only)
  
  # 2. Construct ONLY the filename
  filename <- paste0(prefix, safe_name, ".png")
  
  # 3. Join with the folder path at the final step
  final_path <- file.path(folder, filename)

  p <- ggplot(data, aes(x = display_label, y = mean_value)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
    # Setting base_size here scales everything (axes, titles, etc.)
    theme_bw(base_size = 25) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = title, x = "Individual", y = "Mean Normalization Value") +
    scale_x_discrete(drop = !keep_all_levels)
  
  # Save
  ggsave(final_path, plot = p, width = 16, height = 10, dpi = 300)
  cat("  Saved:", final_path, "\n")
}

cat("\nGenerating plots...\n\n")

# Plot SCG data
cat("Generating SCG plot...\n")
scg_plot_data <- all_means[all_means$type == "SCG", ]
create_barplot(scg_plot_data, 
               title = "Mean SCG Values per Individual", 
               name_only = "scg_barplot",
               folder = opt$output,
               prefix = "Summary_",
                keep_all_levels = TRUE)

# Plot all TEs
cat("Generating All TEs plot...\n")
te_all_plot_data <- all_means[all_means$type == "All_TEs", ]
create_barplot(
    data = te_all_plot_data, 
    title = "Mean TE Values per Individual (All Families)", 
    name_only = "all_TEs_barplot",
    folder = opt$output,
    prefix = "Summary_",
    keep_all_levels = TRUE)

# Plot each TE class (LINE, DNA, LTR, etc.)
cat("Generating TE class plots...\n")
for (cls in te_classes_unique) {
  cat("\n=== Processing TE class:", cls, "===\n")
  te_class_data <- all_means[all_means$type == cls & all_means$category == "TE_class", ]
  
  cat("Filtered data info:\n")
  cat("  Rows:", nrow(te_class_data), "\n")
  cat("  Columns:", paste(colnames(te_class_data), collapse=", "), "\n")
  cat("  First few rows:\n")
  print(head(te_class_data))
  cat("  mean_value class:", class(te_class_data$mean_value), "\n")
  cat("  display_label class:", class(te_class_data$display_label), "\n")
  
  safe_filename <- gsub("[/?*:|\"<>\\\\]", "_", cls)
  create_barplot(
    data = te_class_data, 
    title = paste("Mean TE Values per Individual -", cls, "Class"),
    name_only = paste0(safe_filename, "_barplot"),
    folder = opt$output,
    prefix = "Class_",
    keep_all_levels = TRUE)
}

# Plot each TE family separately
cat("Generating TE family plots...\n")
for (fam in te_families_unique) {
  te_family_data <- all_means[all_means$type == fam & all_means$category == "TE_family", ]
  safe_filename <- gsub("[/?*:|\"<>\\\\]", "_", fam)
  create_barplot(
    data = te_family_data, 
    title = paste("Mean TE Values per Individual -", fam),
    name_only = paste0(safe_filename, "_barplot"),
    folder = opt$output,
    prefix = "Family_",
    keep_all_levels = TRUE)
}

cat("\nGenerating individual TE line plots (Prefix: TE_)...\n")

# Iterate through every row in the original TE data
for (i in 1:nrow(te_data)) {
  row_id <- rownames(te_data)[i]
  
  # Extract values for this specific row and convert to a data frame for ggplot
  # We select all columns except the 'family' and 'class' columns we added earlier
  values <- as.numeric(te_data[i, !colnames(te_data) %in% c("family", "class")])
  
  # Create a data frame that matches the format expected by the function
  line_data <- data.frame(
    individual = colnames(te_data)[!colnames(te_data) %in% c("family", "class")],
    mean_value = values,
    stringsAsFactors = FALSE
  )
  
  # Apply the ordering factors
  line_data$individual <- factor(line_data$individual, levels = individual_order)
  line_data$display_label <- factor(individual_labels[as.character(line_data$individual)], 
                                    levels = individual_labels[individual_order])
  
  # Call the function with the TE_ prefix
  create_barplot(
    data      = line_data,
    title     = paste("Coverage for:", row_id),
    name_only = row_id,
    folder    = opt$output,
    prefix    = "TE_",
    keep_all_levels = TRUE
  )
}

# Create a combined plot with TEs and SCG
cat("Generating comparison plot...\n")
combined_data <- all_means[all_means$type %in% c("SCG", "All_TEs"), ]
combined_data$type <- factor(combined_data$type, levels = c("SCG", "All_TEs"))

p_combined <- ggplot(combined_data, aes(x = display_label, y = mean_value, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  theme_bw(base_size = 25) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(title = "Comparison: SCG vs All TEs",
       x = "Individual",
       y = "Mean Normalization Value") +
  scale_fill_manual(values = c("SCG" = "steelblue", "All_TEs" = "coral"),
                    labels = c("SCG", "All TEs")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(file.path(opt$output, "Summary_scg_vs_TEs_comparison.png"), 
       plot = p_combined, width = 12, height = 6, dpi = 300)
cat("  Saved:", file.path(opt$output, "scg_vs_TEs_comparison.png"), "\n")

# Save summary table
summary_table <- all_means %>%
  select(category, type, individual, display_label, mean_value) %>%
  arrange(category, type, individual)

write.table(summary_table, file.path(opt$output, "summary_means.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved:", file.path(opt$output, "summary_means.tsv"), "\n")

# Print summary statistics
cat("\n=== Summary Statistics ===\n\n")
cat("Number of SCG entries:", nrow(scg_data), "\n")
cat("Number of TE entries:", nrow(te_data), "\n")
cat("TE classes found:", length(te_classes_unique), "\n")
cat("TE classes:", paste(sort(te_classes_unique), collapse=", "), "\n")
cat("TE families found:", length(te_families_unique), "\n")
cat("TE families:", paste(sort(te_families_unique), collapse=", "), "\n")
cat("\nAll plots saved in '", opt$output, "' directory\n", sep="")