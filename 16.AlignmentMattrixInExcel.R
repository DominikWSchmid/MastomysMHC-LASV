library(Biostrings)
library(openxlsx)

# Load the FASTA file as amino acid sequences
alleles <- readAAStringSet("GuineaNigeria_filtered_noGaps_translation.fasta")

# Convert AAStringSet to character vector
seq_names <- names(alleles)
sequences <- as.character(alleles)

# Split sequences into a list of individual characters
split_sequences <- strsplit(sequences, "")

# Reference sequence (first sequence)
ref_sequence <- split_sequences[[1]]

# Create the matrix
result_matrix <- t(sapply(split_sequences, function(seq) {
  ifelse(seq == ref_sequence, ".", seq) # "." for identical positions, else actual value
}))

# Ensure the first row shows the full reference sequence
result_matrix[1, ] <- ref_sequence

# Add position numbers as the header
positions <- 1:ncol(result_matrix)
colnames(result_matrix) <- as.character(positions)

# Add row names for allele names
rownames(result_matrix) <- seq_names

# Save as Excel file
write.xlsx(result_matrix, file = "alleles_matrix_no_gaps.xlsx", rowNames = TRUE)
