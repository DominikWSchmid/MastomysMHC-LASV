setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")

library(dplyr)
library(readr)
# Use stringr for easier regex handling
library(stringr)


rm(list=ls())


c <-read.csv("GuineaNigeria_allelereport_XL.csv", header = TRUE)


c$ID <- as.factor(c$ID)  




c$SUM <- ave(c$COUNT, c$ID, FUN=sum)
c$RPAF <- c$COUNT/c$SUM #relative proportion of reads retrived for one particular allele 
#out of the total amount of reads returned within one given amplicon (relative per amplicon frequency)

c$ID <- as.factor(c$ID)

str(c)


nID<-nlevels(c$ID)


c$ID <- gsub("_rep", "-rep", c$ID)
c$ID <- gsub("_j", "-j", c$ID)
c$ID <- gsub("-rep2", "-j", c$ID)



suffixes <- c("", "-rep", "-j", "-rep-j")




# Step 1: Remove suffixes temporarily to group by base sample ID
# Extract the base sample ID (before suffix)
c$Base_ID <- gsub("-rep|-j|-rep-j", "", c$ID)

# Modify the Base_ID column to capitalize the first three letters, keeping the rest intact
c <- c %>%
  mutate(Base_ID = str_replace(Base_ID, "^([a-z]{3})", toupper))


# Normalize numeration by zero-padding numbers
c$ID <- gsub("(\\D)-(\\d)(?!\\d)", "\\1-0\\2", c$ID, perl = TRUE)


# Extract base IDs without suffixes for grouping
c$Base_ID <- gsub("(-rep|-j|-rep-j)$", "", c$ID)







# Step 2: Summarize total COUNT per replicate for each Base_ID
# This will calculate the total COUNT for each specific ID (with suffixes)
sum_counts <- c %>%
  group_by(Base_ID, ID) %>%
  summarise(Total_COUNT = sum(COUNT, na.rm = TRUE), .groups = 'drop')


# Step 3: Identify the ID with the highest COUNT per Base_ID
# For each base ID, find the replicate with the highest Total_COUNT
# Ensure the data is in a data frame
sum_counts <- as.data.frame(sum_counts)
max_counts <- sum_counts %>%
  group_by(Base_ID) %>%
  top_n(1, Total_COUNT) %>%  # Get the row with the highest Total_COUNT per Base_ID
  ungroup()


# Step 4: Use this to filter the original dataframe 'c'
# Keep only rows corresponding to the IDs with highest counts
c_subset <- c %>%
  filter(ID %in% max_counts$ID) %>%
  select(-Base_ID)  # Optionally remove the Base_ID column

# Result: c_subset now contains only one replicate per sample, the one with the highest sum of COUNT



write_csv(c_subset, "Nigeria_Guinea_highestRep_XL.csv")




########### Examine Data 
# Step 1: Create the Base_ID by removing suffixes
c$Base_ID <- gsub("-rep|-j|-rep-j", "", c$ID)

# Step 2: Summarize total COUNT for each ID within each Base_ID
sum_counts <- c %>%
  group_by(Base_ID, ID) %>%
  summarise(Total_COUNT = sum(COUNT)) %>%
  ungroup()

# Step 3: Identify the highest-count ID for each Base_ID
max_counts <- sum_counts %>%
  group_by(Base_ID) %>%
  filter(Total_COUNT == max(Total_COUNT)) %>%
  slice(1) %>%  # Select the first in case of ties
  rename(Highest_COUNT_ID = ID) %>%
  ungroup()

# Step 4: Gather all replicate IDs except the highest-count one
all_replicates <- sum_counts %>%
  left_join(max_counts, by = "Base_ID") %>%
  filter(ID != Highest_COUNT_ID) %>%
  group_by(Base_ID) %>%
  summarise(All_Other_Replicates = paste(ID, collapse = ", "))

# Step 5: Combine the highest count ID and other replicates into a single table
result_table <- max_counts %>%
  left_join(all_replicates, by = "Base_ID") %>%
  select(Base_ID, Highest_COUNT_ID, All_Other_Replicates)

# View the result
View(result_table)




# Step 1: Extract ALLELE values from both dataframes
alleles_c <- unique(c$ALLELE)         # Unique ALLELEs in the original dataframe
alleles_c_subset <- unique(c_subset$ALLELE)  # Unique ALLELEs in the subset

# Step 2: Identify ALLELEs in 'c' that are not in 'c_subset'
alleles_not_in_subset <- setdiff(alleles_c, alleles_c_subset)

# Convert to a dataframe for better viewing if needed
result_df <- data.frame(Missing_ALLELE = alleles_not_in_subset)

# View the result
print(result_df)

