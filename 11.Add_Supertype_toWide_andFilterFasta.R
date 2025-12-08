setwd("C:/Users/Jan_S/Desktop/2024_11_12_Post-Acacia-Final")


W <- read.csv("Nigeria_Guinea_Wide.csv")

ST <- read.csv("19supertypes.long.csv") 

library(readr)
library(dplyr)
library(Biostrings)
str(W)

str(ST)



# Ensure column names in W match alleleID format in ST
colnames(W) <- gsub("\\.", "_", colnames(W))

# Step 1: Map supertypes to alleles
# Filter columns in W that correspond to ManaMHCI alleles
allele_columns <- grep("^ManaMHCI_", colnames(W), value = TRUE)

# Create a mapping for alleles to supertypes
allele_to_supertype <- setNames(ST$supertype, ST$alleleID)

# Step 2: Add binary columns for each supertype
supertypes <- unique(ST$supertype)  # List of unique supertypes

for (stype in supertypes) {
  # Get alleles corresponding to the current supertype
  alleles_in_supertype <- names(allele_to_supertype[allele_to_supertype == stype])
  
  # Check if the alleles are present in W's columns
  relevant_columns <- intersect(alleles_in_supertype, allele_columns)
  
  # Create a binary column for the supertype
  W[[paste0("Supertype_", stype)]] <- rowSums(W[, relevant_columns, drop = FALSE]) > 0
  W[[paste0("Supertype_", stype)]] <- as.integer(W[[paste0("Supertype_", stype)]])
}

# View the updated dataframe
head(W)
colnames(W)


W$number_alleles <- rowSums(W[, grepl("^ManaMHCI_", colnames(W))] == 1)
W$number_ST <- rowSums(W[, grepl("^Supertype_", colnames(W))] == 1)




W <- W %>% select(1:which(colnames(W) == "Run"),
                  number_alleles, number_ST,
                  everything())

W <- W %>%
  mutate(State = ifelse(Village %in% c("Sokourala", "Sonkonia"), "Faranah", State))

write.csv(W, "Nigeria_Guinea_Wide_ST.csv", row.names = FALSE)


############ Filter 
length(W)

names(W) <- gsub("^ManaMHCI_", "ManaMHC_", names(W))  # Standardize column names

fasta <- readDNAStringSet("GuineaNigeria_filtered_noGaps.fasta")

str(fasta)

alleles_in_fasta <- names(fasta)



# Step 2: Identify ManaMHCI columns that match alleles in the FASTA
mana_mhci_cols <- grep("^ManaMHCI", names(W), value = TRUE)
mana_mhci_cols_to_keep <- intersect(mana_mhci_cols, alleles_in_fasta)

# Step 3: Subset the dataframe
W_filtered <- W %>%
  select(-one_of(setdiff(mana_mhci_cols, mana_mhci_cols_to_keep)))




write.csv(W_filtered, "Nigeria_Guinea_Wide_ST_noGaps.csv", row.names = FALSE)

