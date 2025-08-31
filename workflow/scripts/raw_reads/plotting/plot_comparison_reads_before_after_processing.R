# workflow/scripts/plot_read_counts.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

plot_read_counts <- function(df, species, target_file) {
  file_path <- file.path(target_file)

  df_sum <- df %>%
    summarise(
      raw_count = sum(raw_count, na.rm = TRUE),
      adapter_removed_count = sum(adapter_removed_count, na.rm = TRUE),
      quality_filtered_count = sum(quality_filtered_count, na.rm = TRUE)
    )

  df_long <- df_sum %>%
    pivot_longer(cols = everything(),
                 names_to = "read_type", values_to = "count") %>%
    mutate(read_type = factor(read_type,
                              levels = c("raw_count",
                                         "adapter_removed_count",
                                         "quality_filtered_count")))

  p <- ggplot(df_long, aes(x = species, y = count, fill = read_type)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = c("raw_count" = "#1f77b4",
                                 "adapter_removed_count" = "#ff7f0e",
                                 "quality_filtered_count" = "#2ca02c")) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Species", y = "Read Count (summed)", fill = "Read Type") +
    theme_bw() +
    theme(panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1))

  ggsave(file_path, plot = p, width = 8, height = 5, dpi = 300)
}

plot_compare_reads_before_after_processing <- function(species, source_file, target_file) {
  df <- read.csv(source_file, header = TRUE)
  plot_read_counts(df, species, target_file)
}

# Snakemake integration
plot_compare_reads_before_after_processing(
  snakemake@wildcards$species,
  snakemake@input[[1]],
  snakemake@output[[1]]
)
