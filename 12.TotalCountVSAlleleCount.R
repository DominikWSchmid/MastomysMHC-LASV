setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")

library(readr)
library(dplyr)


XL <- read_csv("Nigeria_Guinea_XL_filtered.csv")

XL$Sample <- gsub("-rep-j|-rep|-j", "", XL$ID)

# Summarize the COUNT and unique ALLELE counts for each ID
summary_df <- XL %>%
  group_by(ID, Sample, Run) %>%
  summarize(
    Total_COUNT = sum(COUNT),
    Allele_Count = n_distinct(ALLELE)
  ) %>%
  ungroup()

# Create a list with each Sample and required information
sample_info_list <- summary_df %>%
  group_by(Sample) %>%
  summarise(
    Highest_Total_COUNT_ID = ID[which.max(Total_COUNT)],          # ID with highest Total_COUNT
    Highest_Allele_Count_ID = ID[which.max(Allele_Count)],        # ID with highest Allele_Count
    All_IDs = list(unique(ID))                                    # List of all IDs in the Sample
  ) %>%
  ungroup() %>%
  group_split(Sample)  # Split into a list by Sample

# Name each element of the list with the Sample name
names(sample_info_list) <- unique(summary_df$Sample)

# Example: Display the first element in the list to verify
sample_info_list[[1]]


# Create the dataframe with an additional column to indicate if the IDs match
sample_info_df <- summary_df %>%
  group_by(Sample) %>%
  summarise(
    Highest_Total_COUNT_ID = ID[which.max(Total_COUNT)],          # ID with the highest Total_COUNT
    Highest_Allele_Count_ID = ID[which.max(Allele_Count)],        # ID with the highest Allele_Count
    All_IDs = paste(unique(ID), collapse = ", "),                 # Concatenate all IDs in the Sample
    Match_Flag = ifelse(ID[which.max(Total_COUNT)] == ID[which.max(Allele_Count)], 1, 0) # 1 if IDs match, 0 if they don't
  ) %>%
  ungroup()

# Display the resulting dataframe
sample_info_df





match_counts <- sample_info_df %>%
  count(Match_Flag)

# Display the result
match_counts







# Add columns for the Total_COUNT of both the ID with the highest Total_COUNT and the ID with the highest Allele_Count
sample_info_df <- summary_df %>%
  group_by(Sample) %>%
  summarise(
    Highest_Total_COUNT_ID = ID[which.max(Total_COUNT)],                # ID with the highest Total_COUNT
    Highest_Total_COUNT_Value = max(Total_COUNT),                       # Total_COUNT of the highest Total_COUNT ID
    Highest_Allele_Count_ID = ID[which.max(Allele_Count)],              # ID with the highest Allele_Count
    Highest_Allele_COUNT_Value = Total_COUNT[which.max(Allele_Count)],  # Total_COUNT of the highest Allele_Count ID
    All_IDs = paste(unique(ID), collapse = ", "),                       # Concatenate all IDs in the Sample
    Match_Flag = ifelse(
      ID[which.max(Total_COUNT)] == ID[which.max(Allele_Count)], 1, 0   # 1 if IDs match, 0 if they don't
    )
  ) %>%
  ungroup()

# Display the updated dataframe
sample_info_df