### Takes in BLAT output and parses it, leaving the "top hit" for each query based on e-value. 
### assumes BLAT command used "-out=blast9"

# Read command line arguments
library(phylotools)
args <- commandArgs(trailingOnly = TRUE)

# Check if there are enough arguments
if (length(args) < 3) {
    cat("Usage: Rscript your_r_script.R input_file sample_name output_directory\n")
    quit(status = 1)
}

# Get input file name, sample name, and output directory from command line arguments
input_file <- args[1]
sample <- args[2]
output_dir <- args[3]

# Construct output file paths based on sample name and output directory
outfile <- file.path(output_dir, paste0(sample, "_IGK.vs.IMGT.txt"))
outfile2 <- file.path(output_dir, paste0(sample, "_IGK.vs.IMGT_consolidated.txt"))

# Read the input BLAST result file
x <- read.table(input_file, as.is = TRUE, sep = "\t")

# Extract the "top hit" for each query based on e-value
out <- tapply(1:nrow(x), x[, 1], function(t) {
  m <- min(x[t, 11])
  y <- x[t, ]
  y[y[, 11] == m, ]
})
kk <- do.call(rbind, out)

#colnames(kk) <- c("Query_id", "Subject_id", "alignment_length", "perc_identity", "mismatches", "gap_openings", "q_star#t", "q_end", "s_start", "s_end", "e.value", "bit_score", "Query_seq", "Subject_seq")

colnames(kk) <- c("Query_id", "Subject_id", "perc_identity", "alignment_length", "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", "e.value", "bit_score")

# Add a new column for the sample name
kk$Sample <- sample

fasta_path <- file.path("/home/egenge01/projects/IGK/run_all_samples/alleles_from_contigs", 
                        sample, 				
                        paste0(sample, "_allele_seqs.fasta"))
# Read the FASTA file
sequence_df <- read.fasta(file = fasta_path, clean_name = FALSE)
# Rename the columns for clarity
colnames(sequence_df) <- c("Query_id", "Query_seq")
# Define the path to the IGK alleles FASTA file
igk_fasta_path <- "/home/egenge01/anaconda3/envs/IGv2/lib/python2.7/site-packages/IGenotyper-1.1-py2.7.egg/IGenotyper/data/alleles.fasta"
# Read the IGK alleles FASTA file
igk_sequence_df <- read.fasta(file = igk_fasta_path, clean_name = FALSE)
# Rename the columns for clarity
colnames(igk_sequence_df) <- c("Subject_id", "Subject_seq")
# Merge the data frames based on matching values of "Query_id" and "Subject_id"
kk <- merge(kk, sequence_df, by = "Query_id", all.x = TRUE)
kk <- merge(kk, igk_sequence_df, by = "Subject_id", all.x = TRUE)

# Extract haplotype from Query_id
kk$sample_haplotype <- ifelse(grepl("hap1", kk$Query_id), paste0(sample, "_1"), paste0(sample, "_2"))

# Extract the allele number from Subject_id
kk$allele_number <- gsub(".*allele=(\\d+).*", "\\1", kk$Subject_id)

# Extract the gene name from Query_id
kk$gene_name <- sub(".*_([^_]+)$", "\\1", kk$Query_id)

# Extract chrom, start, and end from Query_id
chrom_start_end <- gsub(".*chr2:(\\d+)-(\\d+)_hap.*", "chr2 \\1 \\2", kk$Query_id)
chrom_start_end <- strsplit(chrom_start_end, " ")
chrom_start_end <- matrix(unlist(chrom_start_end), ncol = 3, byrow = TRUE, dimnames = list(NULL, c("chrom", "start", "end")))
kk <- cbind(kk, chrom_start_end)

# Write the "top hit" data to the first output file
write.table(kk, outfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Consolidate results for each unique query
uni <- unique(kk$Query_id)
All <- NULL
for (i in uni) {
  sub <- subset(kk, kk$Query_id == i)
  se <- as.vector(sub$Subject_id)
  seqs <- paste(se, collapse = "; ")
  cr <- data.frame(Query_id = i)
  cr$Subject_id <- ifelse(sub$perc_identity[1] == 100, seqs, "NOVEL")
  cr$Sample <- sample
  cr$sample_haplotype <- ifelse(grepl("hap1", i), paste0(sample, "_1"), paste0(sample, "_2"))
  
  # Extract the allele number from Subject_id
  cr$allele_number <- ifelse(sub$perc_identity[1] == 100, gsub(".*allele=(\\d+).*", "\\1", sub$Subject_id[1]), "NOVEL")
  
  cr$gene_name <- sub(".*_([^_]+)$", "\\1", i)
  cr$chrom <- sub(".*chr(\\d+):.*", "chr\\1", i)
  cr$start <- sub(".*chr\\d+:(\\d+)-(\\d+)_.*", "\\1", i)
  cr$end <- sub(".*chr\\d+:\\d+-(\\d+)_.*", "\\1", i)

  cr$perc_identity <- sub(".*\\t([^\\t]+)$", "\\1", sub$perc_identity[1])
  cr$alignment_length <- sub(".*\\t([^\\t]+)$", "\\1", sub$alignment_length[1])
  cr$mismatches <- sub(".*\\t([^\\t]+)$", "\\1", sub$mismatches[1])
  cr$q_start <- sub(".*\\t([^\\t]+)$", "\\1", sub$q_start[1])
  cr$q_end <- sub(".*\\t([^\\t]+)$", "\\1", sub$q_end[1])
  cr$s_start <- sub(".*\\t([^\\t]+)$", "\\1", sub$s_start[1])
  cr$s_end <- sub(".*\\t([^\\t]+)$", "\\1", sub$s_end[1])
  
  # Retain the Query_seq and Subject_seq values from the input file
  cr$Query_seq <- sub(".*\\t([^\\t]+)$", "\\1", sub$Query_seq[1])
  cr$Subject_seq <- sub(".*\\t([^\\t]+)$", "\\1", sub$Subject_seq[1])
  
  All <- rbind(All, cr)
}

# Write the consolidated results to the second output file
write.table(All, outfile2, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)