suppressPackageStartupMessages(library(shazam))

args <- commandArgs(trailingOnly = TRUE)
inputFileName <- args[1]  # Expecting the full path to the .tsv file as input


# Read the input Change-O file
changeoTable <- read.table(inputFileName, sep="\t", header=TRUE)

# Calculate the threshold using Hamming distance
dist_ham <- distToNearest(changeoTable, sequenceColumn="junction", vCallColumn="v_call", jCallColumn="j_call", model="ham", normalize="len", nproc=5)
output <- findThreshold(dist_ham$dist_nearest, method="density")
threshold <- output@threshold
print(threshold)

# Run DefineClones.py with the calculated threshold
system(paste0("DefineClones.py -d ", inputFileName," --act set --nproc 6 --model ham --norm len --dist ", threshold), intern = FALSE, ignore.stdout = FALSE)
