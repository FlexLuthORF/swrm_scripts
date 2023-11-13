library(ggplot2)
library(dplyr)


# Define command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) != 2) {
  cat("Usage: Rscript script.R sample_id root_dir\n")
  quit(status = 1)
}

# Extract sample_id and root_dir from command line arguments
sample_id <- args[1]
root_dir <- args[2]

# Define input file path based on arguments
input_file_path <- file.path(root_dir, "presto", sample_id, "S5_db-pass_clone-pass_germ-pass.tsv")

# Load data
data <- read.csv(input_file_path, sep="\t")

# Basic data processing
data_processed <- data %>%
  filter(!is.na(c_call)) %>%
  group_by(c_call) %>%
  summarise(
    Count = n(),
    WeightedCountConsensus = sum(consensus_count),
    WeightedCountDuplicate = sum(duplicate_count)
  )

# Print the tables
cat("Table 1: Isotype Distribution by Unique Read Count\n")
print(data_processed)

# Save the tables as CSV files with sample_id prefix
table_output_path <- file.path(root_dir, paste0(sample_id, "_isotype_distribution.csv"))
write.csv(data_processed, table_output_path, row.names = FALSE)

# Calculate totals
total_unique_count <- sum(data_processed$Count)
total_weighted_consensus <- sum(data_processed$WeightedCountConsensus)
total_weighted_duplicate <- sum(data_processed$WeightedCountDuplicate)

cat("\nTotal Unique Read Count: ", total_unique_count, "\n")
cat("Total Weighted by Consensus Count: ", total_weighted_consensus, "\n")
cat("Total Weighted by Duplicate Count: ", total_weighted_duplicate, "\n")

# Save the totals as a CSV file with sample_id prefix
totals_output_path <- file.path(root_dir, paste0(sample_id, "_totals.csv"))
totals <- data.frame(
  Metric = c("Total Unique Read Count", "Total Weighted by Consensus Count", "Total Weighted by Duplicate Count"),
  Value = c(total_unique_count, total_weighted_consensus, total_weighted_duplicate)
)
write.csv(totals, totals_output_path, row.names = FALSE)

# Function to create a pie chart
create_pie_chart <- function(data, weight_column, title) {
  ggplot(data, aes(x="", y=!!as.name(weight_column), fill=c_call)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    theme_bw() +  # Change background to white
    labs(fill="Isotype", title=title) +
    geom_text(aes(label=sprintf("%s\n%.1f%%\n(%d)", c_call, (WeightedCountConsensus/total_weighted_consensus)*100, WeightedCountConsensus)),
              position = position_stack(vjust=0.5)) +  # Add labels with percentage and count
    guides(fill=guide_legend(title="Isotype"))  # Add legend label
}

# Create pie charts
pie1 <- create_pie_chart(data_processed, "Count", "Isotype Distribution by Unique Read Count")
pie2 <- create_pie_chart(data_processed, "WeightedCountConsensus", "Isotype Distribution Weighted by Consensus Count")
pie3 <- create_pie_chart(data_processed, "WeightedCountDuplicate", "Isotype Distribution Weighted by Duplicate Count")

# Save the plots with sample_id prefix
ggsave(file.path(root_dir, "presto/analysis", paste0(sample_id, "_pie_unique.png")), pie1)
ggsave(file.path(root_dir, "presto/analysis", paste0(sample_id, "_pie_consensuscount.png")), pie2)
ggsave(file.path(root_dir, "presto/analysis", paste0(sample_id, "_pie_duplicates.png")), pie3)
