# Load necessary libraries
library(airr)
library(shazam)
library(dplyr)
library(alakazam)

# Read the rearrangement data
airr <- read_rearrangement("changeo/S5_filter-pass_germ-pass.tsv")

# Step 8: Identify and Remove Chimeric Sequences
is_chimeric <- slideWindowDb(
  airr,
  sequenceColumn = "sequence_alignment",
  germlineColumn = "germline_alignment_d_mask",
  mutThresh=6,
  windowSize=10
)

# Print the table of chimeric sequences
print(table(is_chimeric))

# Remove chimeric sequences
airr <- airr[!is_chimeric,]

# Step 9: Collapse Duplicates
num_fields <- c("consensus_count", "duplicate_count")
collapse_groups <- c("v_gene", "j_gene", "junction_length", "c_call", "productive")

airr <- airr %>%
  mutate(v_gene=getGene(v_call), j_gene=getGene(j_call)) %>%
  group_by(.dots=collapse_groups) %>%
  do(collapseDuplicates(.,
    id = "sequence_id",
    seq = "sequence_alignment",
    text_fields = NULL,
    num_fields = num_fields,
    seq_fields = NULL,
    add_count = TRUE,
    ignore = c("N", "-", ".", "?"),
    sep = ",",
    dry = FALSE,
    verbose = FALSE
  )) %>%
  ungroup() %>%
  select(-v_gene, -j_gene)

# Optionally, save the final result to a new file
write_rearrangement(airr, file="changeo/S5_filter-pass_germ-pass_collapsed.tsv")

