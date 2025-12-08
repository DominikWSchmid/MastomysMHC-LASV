
setwd("C:/Users/Jan_S/Desktop/2024_11_12_Post-Acacia-Final")


setwd("C:/Users/jansa/Desktop/Ayo/2024_11_12_Post-Acacia-Final")

library(tidyverse)

library(dplyr)
library(readxl)
library(rlang)

library(ggplot2)

# Started out in the lab with the following 

NigeriaSamp <- scan("SamplesNigeria2024.txt", what = "", sep = ",")

Nigeriasamples <- data.frame(ID = NigeriaSamp) 



GuineaSamp <- scan("SamplesGuinea2024.txt", what = "", sep = ",")

Guineasamples <- data.frame(ID = GuineaSamp) 

# Step 1: Ensure first three characters (prefix) are capital letters
Guineasamples$ID <- sub("^(\\w)(\\w+)", "\\U\\1\\L\\2", Guineasamples$ID, perl = TRUE)

# Step 2: Replace the blank space with an underscore
Guineasamples$ID <- gsub(" ", "-", Guineasamples$ID)

# Step 3: Ensure that the numbers are three digits
Guineasamples$ID <- gsub("(\\D)(\\d{1})(\\D|$)", "\\100\\2\\3", Guineasamples$ID)
Guineasamples$ID <- gsub("(\\D)(\\d{2})(\\D|$)", "\\10\\2\\3", Guineasamples$ID)



str(Guineasamples)
str(Nigeriasamples)

### metadata 

os <- read_excel("MHC_OneSheet.xlsx", na = "NA")
os$Village[os$Village == "Eguare-Egoro"] <- "Ekpoma"
names(os)[names(os) == "Sample_No."] <- "ID"
os <- os


gm <- read_excel("MHC_Guinea.xlsx", na = "NA")
colnames(gm)[colnames(gm) == "Nviro"] <- "ID"

str(os)
str(gm)



# Step 1: Ensure first three characters (prefix) are capital letters
gm$ID <- sub("^(\\w)(\\w+)", "\\U\\1\\L\\2", gm$ID, perl = TRUE)
# Step 2: Replace the blank space with an underscore
gm$ID <- gsub(" ", "-", gm$ID)
# Step 3: Ensure that the numbers are three digits
gm$ID <- gsub("(\\D)(\\d{1})(\\D|$)", "\\100\\2\\3", gm$ID)
gm$ID <- gsub("(\\D)(\\d{2})(\\D|$)", "\\10\\2\\3", gm$ID)



gm$Village <- NA
gm <- gm %>%
  mutate(Village = case_when(
    grepl("^Dam", ID) ~ "Damania",
    grepl("^Tam", ID) ~ "Tambaya",
    grepl("^Son", ID) ~ "Sonkonia",
    grepl("^Sok", ID) ~ "Sokourala",
    TRUE ~ Village  # Keep the original Village value if no match
  ))
# What has gone through Acacia, ergo has been sequenced


PostseqAyo <- read.csv("Ayoold/Ayoold_pipelinereport.csv")

PostseqG <- read.csv("Guinea_pipelinereport.csv")                  

PostseqN <- read.csv("Nigeria_pipelinereport.csv")        



# Step 1: Ensure first three characters (prefix) are capital letters
PostseqG$ID <- sub("^(\\w)(\\w+)", "\\U\\1\\L\\2", PostseqG$ID, perl = TRUE)
# Step 2: Replace the blank space with an underscore
PostseqG$ID <- gsub(" ", "-", PostseqG$ID)
# Step 3: Ensure that the numbers are three digits
PostseqG$ID <- gsub("(\\D)(\\d{1})(\\D|$)", "\\100\\2\\3", PostseqG$ID)
PostseqG$ID <- gsub("(\\D)(\\d{2})(\\D|$)", "\\10\\2\\3", PostseqG$ID)



# Filtered for PAML

Supertyped <- read.csv("Nigeria_Guinea_Wide_ST.csv")


# Create Sample-ID
Supertyped$Sample_ID <- Supertyped$ID 


PostseqAyo$Sample_ID <- gsub("(-rep)", "", PostseqAyo$ID)

PostseqG$Sample_ID <- gsub("(-j|-rep|-rep-j)$", "", PostseqG$ID)


PostseqN$ID <- gsub("^NG-", "", PostseqN$ID)
PostseqN$Sample_ID <- gsub("(-j|-rep|-rep-j)$", "", PostseqN$ID)


# Step 1: Add Country Information
os <- os %>% mutate(Country = "Nigeria", ID = as.character(ID))
gm <- gm %>% mutate(Country = "Guinea", ID = as.character(ID))

# Step 2: Combine Metadata from Nigeria and Guinea
metadata <- bind_rows(
  os %>% select(ID, Country, Village),
  gm %>% select(ID, Country, Village)
)

# Step 3: Combine All Sample Datasets
combined_samples <- bind_rows(
  PostseqAyo %>% select(ID, Sample_ID),
  PostseqG %>% select(ID, Sample_ID),
  PostseqN %>% select(ID, Sample_ID)
)
 


# Combine all datasets into one
PostseqCombined <- rbind(
  PostseqAyo[, c("ID", "Sample_ID")],
  PostseqG[, c("ID", "Sample_ID")],
  PostseqN[, c("ID", "Sample_ID")]
)


# Step 4: Merge Combined Samples with Metadata
combined_samples_with_metadata <- combined_samples %>%
  left_join(metadata, by = c("Sample_ID" = "ID"))

# Correct way to count occurrences of each Sample_ID
combined_samples_with_metadata %>%
  group_by(Sample_ID) %>%
  tally() %>%
  filter(n > 1)


# Check the count of each Sample_ID
combined_samples_with_metadata %>%
  group_by(Sample_ID) %>%
  summarise(count = n())

combined_samples_with_metadata %>%
  group_by(Sample_ID) %>%
  summarise(count = dplyr::n(), .groups = "drop")


combined_metadata <- combined_samples_with_metadata %>%
  group_by(Sample_ID) %>%
  summarise(
    ID = paste(unique(ID), collapse = ", "),  # Combine unique IDs as a string
    Country = first(Country),  # Consistent country for a given Sample_ID
    Village = first(Village),  # Consistent village for a given Sample_ID
    Replicate_Count = n_distinct(ID)  # Count unique IDs per Sample_ID
  ) %>%
  ungroup()


# Step 5: Group by Sample_ID and calculate correct Replicate_Count based on unique ID values
combined_metadata <- combined_samples_with_metadata %>%
  group_by(Sample_ID) %>%
  summarise(
    ID = paste(unique(ID), collapse = ", "),  # Combine unique IDs as a string
    # Processed_By = paste(unique(Processed_By[Processed_By != ""]), collapse = ", "),  # Handle empty Processed_By
    Country = first(Country),  # Consistent country for a given Sample_ID
    Village = first(Village),  # Consistent village for a given Sample_ID
    Replicate_Count = n_distinct(ID)  # Count unique IDs per Sample_ID
  ) %>%
  ungroup()


# Step 6: Add Status Flags and Ensure Consistency
combined_metadata <- combined_metadata %>%
  mutate(
    # Processed = Sample_ID %in% PostseqCombined$Sample_ID,
    Sequenced = Sample_ID %in% PostseqCombined$Sample_ID[PostseqCombined$Processed_By %in% c("Nig", "Gui")],
    Supertyped = Sample_ID %in% Supertyped$Sample_ID,
    Sequenced = ifelse(Supertyped & !Sequenced, TRUE, Sequenced)  # Ensure consistency
  )



# Check for remaining inconsistencies
inconsistent_samples <- combined_metadata %>%
  filter(Supertyped & !Sequenced)
print(inconsistent_samples)




# Update Replicate_Count to count the number of unique IDs for each Sample_ID
combined_metadata <- combined_metadata %>%
  mutate(
    Replicate_Count = sapply(strsplit(ID, ",\\s*"), length)  # Split IDs by ', ' and count
  )


# Check if the `Replicate_Count` column was added correctly
head(combined_metadata)



str(combined_samples_with_metadata)
str(combined_metadata)

table(duplicated(combined_metadata$Sample_ID))  # This should ideally be 0 for each Sample_ID
head(combined_metadata %>% group_by(Sample_ID) %>% summarise(Replicate_Count = n()))

############################# report  new #############


# Create summary table grouped by Country
summary_table <- combined_metadata %>%
  group_by(Village) %>%
  summarise(
    Country = first(Country), 
    Total_Sample_IDs = n(),  # Total number of Sample_IDs
    With_Replicates = sum(Replicate_Count >= 2),  # Count of Sample_IDs with at least 2 replicates
    Sequenced_Count = sum(Sequenced),  # Count of Sample_IDs marked as Sequenced
    Supertyped_Count = sum(Supertyped)  # Count of Sample_IDs marked as Supertyped
  )

# Display the summary table
summary_table


write.csv(summary_table, "summary_table.csv", row.names = FALSE)


