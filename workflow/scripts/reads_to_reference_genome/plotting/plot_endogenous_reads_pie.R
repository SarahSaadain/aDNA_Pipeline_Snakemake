library(dplyr)
library(ggplot2)

plot_endogenous_reads_pie <- function(source_file, output_file) {
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  df <- read.table(source_file, sep =",", header = TRUE) %>%
    mutate(
      Non_Endogenous = TotalReads - MappedReads,
      Proportion = MappedReads / TotalReads,
      percent_endogenous = Proportion * 100, 
      percent_non_endogenous = (Non_Endogenous / TotalReads) * 100
    )
  
  df$Individual <- sub("\\..*", "", df$Filename)
  
  # Generate all individual plots into one file (multi-page PDF)
  pdf(output_file, width = 6, height = 6)
  for(i in 1:nrow(df)) {
    row_data <- df[i, ]
    
    row_long <- data.frame(
      read_type = c("Endogenous", "Non-Endogenous"),
      count = c(row_data$MappedReads, row_data$Non_Endogenous),
      percent = c(row_data$percent_endogenous, row_data$percent_non_endogenous)
    )
    
    p <- ggplot(row_long, aes(x = "", y = count, fill = read_type)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      geom_text(aes(label = paste0(round(percent, 1), "%")), 
                position = position_stack(vjust = 0.5), size = 5, color = "white") +
      scale_fill_manual(values = c("Endogenous" = "#209557", "Non-Endogenous" = "#1f5bb4")) +
      labs(x = NULL, y = NULL, fill = "Read Type", title = paste("Endogenous vs Non-Endogenous Reads:", row_data$Individual)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            legend.position = "bottom",
            panel.border = element_blank(),
            plot.border = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
    
    print(p)
  }
  dev.off()
}

if (exists("snakemake")) {
  source_file <- snakemake@input[[1]]
  output_file <- snakemake@output[[1]]
  plot_endogenous_reads_pie(source_file, output_file)
}
