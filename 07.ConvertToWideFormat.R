setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")


# Load necessary library

library(tidyverse)
library(dplyr)
library(data.table)



# Lese die Daten ein
XL <- read.csv("GN_noRep_XL_filtered.csv")

# Überprüfe die Struktur der Daten


data<-data.table(XL)



wide_data <- data %>%
  # First convert the data to wide format
  dcast(ID + Run ~ ALLELE, value.var = "ALLELE", fun.aggregate = length) %>%
  # Convert counts to 0/1 for all columns except ID and Run
  mutate(across(-c(ID, Run), ~ ifelse(. > 0, 1, 0))) %>%
  # Filter out rows where ID is NA (if there are any such rows)
  filter(!is.na(ID))


# Optionally, save to a new CSV file
write.csv(wide_data, "Nigeria_Guinea_wide_allele.csv", row.names = FALSE)
