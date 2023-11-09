# Load necessary libraries
library(airr)
library(alakazam)
library(stringi)
library(dplyr)

# Read the rearrangement data
airr <- read_rearrangement("changeo/S5_db-pass_clone-pass.tsv")

# Step 1: Identify Short Sequences (Minimum length 200 nt)
long_seq <- stri_count(airr[['sequence_alignment']], regex="[^-.N]") >= 200

# Step 2: Identify Reads with Coherent Gene, Primer, and Isotype Calls
same_locus <- getLocus(airr[['v_call']]) == airr[['locus']] & 
              getLocus(airr[['c_call']]) == airr[['locus']]

# Step 3: Identify Reads with an Acceptable Number of Ambiguous Nucleotides (Max 10% N)
num_n <- stri_count(airr[['sequence_alignment']], fixed="N")
len <- stri_count(airr[['sequence_alignment']], regex="[^-.]")
low_n <- num_n / len <= 0.10

# Step 4: Identify Productive Sequences
prod <- airr[['productive']]

# Step 5: Identify Sequences with Junction Length Multiple of Three
m3 <- airr[['junction_length']] %% 3 == 0

# Step 6: Filter and Save
filter_pass <- long_seq & same_locus & low_n & prod & m3
write_rearrangement(airr[filter_pass, ], file="changeo/S5_db-pass_clone-pass_filter-pass.tsv")

