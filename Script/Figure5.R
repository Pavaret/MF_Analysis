# Install required packages if not already installed
if_require <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE, repos = "https://cran.rstudio.com/")
    library(package, character.only = TRUE)
  }
}

required_packages <- c("tidyverse", "ggprism", "ggtext", "ggrepel", "ggthemes", "RColorBrewer", "viridis", "forcats", "readxl", "argparse")
for (pkg in required_packages) {
  if_require(pkg)
}


# Load libraries
library(tidyverse)
library(ggprism)
library(ggtext)
library(ggrepel)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(forcats)
library(readxl)
library(argparse)

# Set up argument parser
parser <- ArgumentParser(description = "Functional Enrichment Plot Generator")
parser$add_argument("--data_dir", type = "character", help = "Directory containing the input data files", required = TRUE)
parser$add_argument("--output_dir", type = "character", help = "Directory to save the output figures", required = TRUE)
args <- parser$parse_args()

# Define a function to perform enrichment analysis and generate plots
create_enrichment_plot <- function(input_file, output_file, title, x_limits, fill_colors) {
  # Load the data
  enrichment_data <- read_xlsx(input_file) %>% 
    mutate(column_0 = str_replace(column_0, "GO:", ""),
           FDR = p.adjust(column_3, method = "BH")) %>%  
    filter(FDR <= 0.05) %>% 
    mutate(Log_FDR = -log(FDR))
  
  # Ordering the GO Modules
  plot_df <- enrichment_data %>%
    group_by(column_0) %>% 
    slice_max(order_by = Log_FDR, n = 15) %>% 
    ungroup() %>% 
    mutate(
      Sub_Cat = fct_relevel(column_0, c("BP", "CC", "MF", "KEGG")),  # Relevel based on given categories
      column_2 = factor(column_2, levels = unique(column_2))  # Reorder Pathway based on arranged order
    ) %>% 
    filter(column_0 %in% c("BP", "CC", "MF", "KEGG"))
  
  # Generate the plot
  plot <- plot_df %>% 
    ggplot(aes(x = Log_FDR, y = column_2, fill = Sub_Cat)) +
    geom_bar(stat = "identity", position = "dodge") +
    guides(fill = guide_legend(ncol = 1, override.aes = list(size = 4))) +
    scale_x_continuous(guide = "prism_minor", expand = c(0.005, 0), breaks = seq(0, x_limits, 5), limits = c(0, x_limits)) +
    labs(x = "-Log(FDR)", y = " ", fill = "GO Term", title = title) +
    scale_fill_manual(values = fill_colors) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 2, face = "bold"),
      title = element_markdown(size = 14, family = "Arial"),
      axis.title = element_markdown(size = 14, family = "Arial", color = "Black"),
      axis.ticks.y = element_line(size = 0),
      axis.text.y = element_markdown(size = 12, family = "Arial", color = "Black"),
      axis.text.x = element_markdown(size = 12, family = "Arial", color = "Black"),
      legend.text = element_markdown(size = 12, family = "Arial"),
      legend.title = element_markdown(size = 14, face = "bold", family = "Arial"),
      legend.key.size = unit(1.2, "line")
    )
  
  # Save the plot
  jpeg(filename = output_file, width = 12, height = 7.5, units = "in", res = 600)
  print(plot)
  dev.off()
}

# Define parameters for Deletion and Insertion plots
fill_colors_del <- c("BP" = "#80558C", "CC" = "#6E85B7", "MF" = "#90C8AC")
fill_colors_ins <- c("BP" = "#80558C", "CC" = "#6E85B7", "MF" = "#90C8AC", "KEGG" = "#D77FA1")

# Create plots for Deletion and Insertion
create_enrichment_plot(
  input_file = file.path(args$data_dir, 'ALL_GPROFILER_DEL_>1K.xlsx'),
  output_file = file.path(args$output_dir, 'Functional_Enrichment_Analysis_of_Deletion_(>1kb).jpg'),
  title = 'Functional Enrichment Analysis of Deletion (>1kb)',
  x_limits = 23,
  fill_colors = fill_colors_del
)

create_enrichment_plot(
  input_file = file.path(args$data_dir, 'ALL_GPROFILER_INS_>1K.xlsx'),
  output_file = file.path(args$output_dir, 'Functional_Enrichment_Analysis_of_Insertion_(>1kb).jpg'),
  title = 'Functional Enrichment Analysis of Insertion (>1kb)',
  x_limits = 43,
  fill_colors = fill_colors_ins
)
