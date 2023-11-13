library(alakazam)
library(airr)
library(dplyr)
library(ggplot2)
library(scales)
color_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

# Function to analyze a single sample
analyze_sample <- function(sample_id, root_dir) {
  file_path <- paste0(root_dir, "/presto/", sample_id, "/S5_db-pass_clone-pass_germ-pass.tsv")
  airr <- read_rearrangement(file_path)

  # Perform your existing analysis here...
  v_usage_isotype <- countGenes(airr, "v_call", group="c_call", fill=TRUE)
  most_used_v <- v_usage_isotype %>%
    filter(c_call != "IGHE") %>%
    group_by(c_call) %>%
    slice_max(., seq_freq, n=1)

  # Plotting
  gene_usage_plot <- ggplot(v_usage_isotype %>%
    filter(c_call != "IGHE" & gene %in% most_used_v[['gene']]),
    aes(x=gene, y=seq_freq, color=c_call)) +
    scale_color_manual(values=color_palette) +
    scale_y_continuous(labels=percent) +
    geom_point(size=2) + theme_bw() +
    xlab("Gene") + ylab("Frequency") +
    guides(color=guide_legend(title="Isotype"))

  # Save plot and table
  plot_path <- paste0(root_dir, "/presto/analysis/", sample_id, "_plot.png")
  table_path <- paste0(root_dir, "/presto/analysis/", sample_id, "_table.tab")

  ggsave(plot_path, gene_usage_plot, width=10, height=6)
  write.table(v_usage_isotype, table_path, sep="\t", row.names=FALSE)
}

# Function to perform combined analysis
# Function to perform combined analysis
analyze_combined <- function(sample_ids, root_dir) {
  # Filter out non-existing files
  valid_sample_ids <- Filter(function(sample_id) {
    file_path <- paste0(root_dir, "/presto/", sample_id, "/S5_db-pass_clone-pass_germ-pass.tsv")
    if (file.exists(file_path)) {
      return(TRUE)
    } else {
      message("Skipping non-existing file for combined analysis: ", file_path)
      return(FALSE)
    }
  }, sample_ids)

  # Check if there are any valid files left
  if (length(valid_sample_ids) == 0) {
    message("No valid files found for combined analysis.")
    return(NULL)
  }

  # Read and combine data from valid files
  combined_airr <- do.call(rbind, lapply(valid_sample_ids, function(sample_id) {
    file_path <- paste0(root_dir, "/presto/", sample_id, "/S5_db-pass_clone-pass_germ-pass.tsv")
    read_rearrangement(file_path)
  }))

  # Perform the combined analysis
  v_usage_isotype_combined <- countGenes(combined_airr, "v_call", group="c_call", fill=TRUE)
  most_used_v_combined <- v_usage_isotype_combined %>%
    filter(c_call != "IGHE") %>%
    group_by(c_call) %>%
    slice_max(., seq_freq, n=1)

  # Plotting for combined data
  gene_usage_plot_combined <- ggplot(v_usage_isotype_combined %>%
    filter(c_call != "IGHE" & gene %in% most_used_v_combined[['gene']]),
    aes(x=gene, y=seq_freq, color=c_call)) +
    scale_color_manual(values=color_palette) +
    scale_y_continuous(labels=percent) +
    geom_point(size=2) + theme_bw() +
    xlab("Gene") + ylab("Frequency") +
    guides(color=guide_legend(title="Isotype"))

  # Save combined plot and table
  plot_path_combined <- paste0(root_dir, "/presto/analysis/combined_plot.png")
  table_path_combined <- paste0(root_dir, "/presto/analysis/combined_table.tab")

  ggsave(plot_path_combined, gene_usage_plot_combined, width=10, height=6)
  write.table(v_usage_isotype_combined, table_path_combined, sep="\t", row.names=FALSE)
}


# Sample IDs and root directory
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  stop("Usage: Rscript script_name.R <path_to_sample_ids_file> <root_dir>")
}

sample_ids_file <- args[1]
root_dir <- args[2]

# Read sample IDs from the file (assuming the first column contains sample IDs)
sample_data <- read.table(sample_ids_file, header = FALSE, sep = "\t")
sample_ids <- sample_data[[1]]

# Analyze each sample
for (sample_id in sample_ids) {
  file_path <- paste0(root_dir, "/presto/", sample_id, "/S5_db-pass_clone-pass_germ-pass.tsv")
  if (file.exists(file_path)) {
    analyze_sample(sample_id, root_dir)
  } else {
    message("Skipping non-existing file: ", file_path)
  }
}

# Perform combined analysis
analyze_combined(sample_ids, root_dir)
