setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")

library(readr)
library(dplyr)
detach("package:plyr", unload=TRUE)###### turn off plyr!!! 
###### turn off plyr!!! 

XL <- read_csv("Nigeria_Guinea_XL_filtered.csv")

# Count the frequency of each ALLELE
allele_counts <- XL %>%
  count(ALLELE)


# Filter the XL dataframe to remove ALLELES that appear only once
XL <- XL %>%
  filter(ALLELE %in% allele_counts$ALLELE[allele_counts$n > 1])



# Remove "-rep-j", "-rep", "-j" suffixes and "NG-" prefix
#XL$ID <- gsub ("^(NG-)","", XL$ID)
XL$Sample <- gsub("^(NG-)|(-rep-j|-rep|-j)", "", XL$ID)


# Summarize the COUNT and unique ALLELE counts for each ID
summary_df <- XL %>%
  group_by(ID, Sample, Run) %>%
  summarize(
    Total_COUNT = sum(COUNT),
    Allele_Count = n_distinct(ALLELE)
  ) %>%
  ungroup()


# Modify to handle ties in Allele_Count by checking if the ID with the highest Total_COUNT is in the set of IDs with the highest Allele_Count
sample_info_df <- summary_df %>%
  group_by(Sample) %>%
  summarise(
    Run = paste(Run)[which.max(Allele_Count)],
    Highest_Total_COUNT_ID = ID[which.max(Total_COUNT)],                # ID with the highest Total_COUNT
    Highest_Total_COUNT_Value = max(Total_COUNT),                       # Total_COUNT of the highest Total_COUNT ID
    Highest_Allele_Count_ID = ID[which.max(Allele_Count)],              # ID with the highest Allele_Count
    Highest_Allele_COUNT_Value = Total_COUNT[which.max(Allele_Count)],  # Total_COUNT of the highest Allele_Count ID
    All_IDs = paste(unique(ID), collapse = ", "),                       # Concatenate all IDs in the Sample
    Num_IDs = n_distinct(ID),                                           # Number of unique IDs in the Sample
    HighestCountHasMostAllele = ifelse(
      Highest_Total_COUNT_ID %in% ID[Allele_Count == max(Allele_Count)], 1, 0
    )  # 1 if the ID with the highest Total_COUNT is among those with the highest Allele_Count
  ) %>%
  ungroup()




# only one ID per sample 

num_ids_1_count <- sample_info_df %>%
  filter(Num_IDs == 1) %>%
  nrow()

# more than one ID per sample
num_ids_more_count <- sample_info_df %>%
  filter(Num_IDs > 1) %>%
  nrow()

# more than one ID, most reads ID has most alleles 
num_ids_highes_count_most_alleles <- sample_info_df %>%
  filter(Num_IDs > 1 & HighestCountHasMostAllele == 1 ) %>%
  nrow()


# more than one ID, most reads ID has NOT most alleles 
num_ids_highes_count_not_alleles <- sample_info_df %>%
  filter(Num_IDs > 1 & HighestCountHasMostAllele == 0 ) %>%
  nrow()




##### 

# Filter out IDs where HighestCountHasMostAllele == 0


GNNG <- read_csv("GuineaNigeria_allelereport_XL.csv")


str(GNNG)
# Get the IDs where HighestCountHasMostAllele == 1
highest_count_ids <- sample_info_df %>%
  filter(HighestCountHasMostAllele == 1) %>%
  pull(Highest_Total_COUNT_ID)

# Filter GNNG to keep only rows with the Highest_Total_COUNT_ID
GNNG_subset <- GNNG %>%
  filter(ID %in% highest_count_ids)



write_csv(GNNG_subset, "GN_noRep_XL.csv")




# Remove "-rep-j", "-rep", "-j" suffixes and "NG-" prefix
#XL$ID <- gsub ("^(NG-)","", XL$ID)
GNNG_subset$Sample <- gsub("^(NG-)|(-rep-j|-rep|-j)", "", GNNG_subset$ID)


# Summarize the COUNT and unique ALLELE counts for each ID
summary_GN <- GNNG_subset %>%
  group_by(ID, Sample, Run) %>%
  summarize(
    Total_COUNT = sum(COUNT),
    Allele_Count = n_distinct(ALLELE)
  ) %>%
  ungroup()


length(unique(GNNG_subset$ID))
