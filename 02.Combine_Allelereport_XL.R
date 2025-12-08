setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")
#combine allelereport_XL

library(readr)
library(dplyr)
library(stringr)

# Function to read and clean CSV files
read_and_clean_csv <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE)
  colnames(df) <- c("ID", "ALLELE", "COUNT", "SEQ")
  return(df)
}

# Read the updated and AYO allele reports
old_XL <- read_and_clean_csv("updated_Nigeria_allelereport_XL.csv")
new_XL <- read_and_clean_csv("updated_Guinea_allelereport_XL.csv")



old_XL$Run <- "NG"
new_XL$Run <- "GN"


# Convert ID columns to character type to avoid type mismatch



old_XL <- old_XL %>%
  mutate(ID = as.character(ID))

new_XL <- new_XL %>%
  mutate(ID = as.character(ID))




# Combine the data frames
combined_sequences <- bind_rows(old_XL, new_XL)


combined_sequences <- combined_sequences %>%
  mutate(ID = str_remove(ID,"^NG-"))

combined_sequences <- combined_sequences %>%
  mutate(ID = str_remove(ID,"^NG-"))

combined_sequences <- combined_sequences[combined_sequences$ID != "PCR-Neg2-j", ]


# Write the combined data frame back to a CSV file
write_csv(combined_sequences, "GuineaNigeria_allelereport_XL.csv")
