# workflow/scripts/plot_read_counts.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)


plot_by_individual <- function(df, species, target_file) {

  file_path <- file.path(target_file)
  
  df_individual <- df %>%
    group_by(individual) %>%
    summarise(
      raw_count_absolute = sum(raw_count),
      adapter_removed_count_absolute = sum(adapter_removed_count),
      quality_filtered_count_absolute = sum(quality_filtered_count),
      .groups = "drop"
    ) %>%

    pivot_longer(cols = c(raw_count_absolute, adapter_removed_count_absolute, quality_filtered_count_absolute),
                 names_to = "read_type", values_to = "count") %>%

    mutate(read_type = factor(read_type,
                              levels = c("raw_count_absolute",
                                         "adapter_removed_count_absolute",
                                         "quality_filtered_count_absolute")))

  p <- ggplot(df_individual, aes(x = individual, y = count, fill = read_type)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = c("raw_count_absolute" = "#1f77b4",
                                 "adapter_removed_count_absolute" = "#ff7f0e",
                                 "quality_filtered_count_absolute" = "#2ca02c")) +
    scale_y_continuous(labels = comma) +
    labs(x = "Individual", y = "Read Count", fill = "Read Type") +
    theme_bw() +
    theme(panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))
          
  ggsave(file_path, plot = p, width = 8, height = 5, dpi = 300)
}

plot_compare_reads_before_after_processing <- function(species, source_file, target_file) {
  df <- read.csv(source_file, header = TRUE)
  plot_by_individual(df, species, target_file)
}

# Snakemake integration
plot_compare_reads_before_after_processing(
  snakemake@wildcards$species,
  snakemake@input[[1]],
  snakemake@output[[1]]
)
