setwd("C:/Users/Jan_S/Desktop/2024_11_12_Post-Acacia-Final")

setwd("C:/Users/jansa/Desktop/Ayo/2024_11_12_Post-Acacia-Final")



library(dplyr)
library(readxl)
library(tidyr)
library(ggplot2)
library (BiocManager)
library (Biostrings)


df <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv")

str(df)



# Summarize data by Country, Village, and LASV_positive
summary_table <- df %>%
  group_by(Country, Village) %>%
  summarise(
    Total_Samples = n(),
    LASV_Positive = sum(LASV_positive == "pos", na.rm = TRUE)
  ) %>%
  arrange(Country, Village)

# Split the summary table by Country into a list of tables
country_tables <- split(summary_table, summary_table$Country)

# View tables for each country
country_tables$Nigeria  # Example for Nigeria
country_tables$Guinea   # Example for Guinea, if applicable


sum(grepl("^ManaMHCI", names(df)))




# Select columns starting with 'ManaMHCI_'
mana_mhci_cols <- df %>%
  select(starts_with("ManaMHCI"))

# Count columns where all values are 0
num_all_zero_cols <- sum(colSums(mana_mhci_cols) == 0)

# Display the result
num_all_zero_cols



# Select columns starting with 'ManaMHCI_'
mana_mhci_cols <- df[, grepl("^ManaMHCI", names(df))]

# Count columns where all values are 0
num_all_zero_cols <- sum(colSums(mana_mhci_cols) == 0)
num_all_zero_cols









fasta <-  readDNAStringSet("GuineaNigeria_filtered_noGaps.fasta")

length(fasta)

NG <-  readDNAStringSet("GuineaNigeria.fasta")
length(NG)


NGf <-  readDNAStringSet("GuineaNigeria_filtered.fasta")
length(NGf)




# Step 1: Identify allele columns (e.g., columns starting with "ManaMHCI")
allele_cols <- grep("^ManaMHCI", names(df), value = TRUE)

# Step 2: Count the number of alleles (non-zero entries) for each ID
allele_counts <- df %>%
  rowwise() %>%
  mutate(Allele_Count = sum(c_across(all_of(allele_cols)) != 0)) %>%
  ungroup()

# Step 3: Compute summary statistics
summary_stats <- allele_counts %>%
  summarise(
    Median = median(Allele_Count),
    Mean = mean(Allele_Count),
    Min = min(Allele_Count),
    Max = max(Allele_Count),
    Q1 = quantile(Allele_Count, 0.25),
    Q3 = quantile(Allele_Count, 0.75)
  )

# Print the results
summary_stats

# Create a histogram
ggplot(allele_counts, aes(x = Allele_Count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "cornflowerblue", alpha = 0.7) +
  labs(
    title = "Distribution of Allele Counts",
    x = "Number of Alleles",
    y = "Frequency"
  ) +
  theme_minimal()


ggplot(allele_counts, aes(x = factor(Allele_Count))) +
  geom_bar(color = "black", fill = "orange", alpha = 0.8) +
  labs(
    title = "Distribution of Allele Counts",
    x = "Number of Alleles",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Tilt x-axis labels if necessary


###############

AllelesNigeria <- subset(allele_counts, allele_counts$Country == "Nigeria")

ggplot(AllelesNigeria, aes(x = Allele_Count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "coral", alpha = 0.7) +
  labs(
    title = "Distribution of Allele Counts in Nigeria",
    x = "Number of Alleles",
    y = "Frequency"
  ) +
  theme_minimal()


AllelesGuinea <- subset(allele_counts, allele_counts$Country == "Guinea")

ggplot(AllelesGuinea, aes(x = Allele_Count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "aquamarine3", alpha = 0.7) +
  labs(
    title = "Distribution of Allele Counts in Guinea",
    x = "Number of Alleles",
    y = "Frequency"
  ) +
  theme_minimal()



################## Villages of INterest 

AllelesInt <- subset(allele_counts, allele_counts$Village %in% c("Damania", "Tambaya", "Sonkonia", 
                                                                 "Sokourala", "Ebudin", "Ekpoma", 
                                                                 "Owo", "Aba gboro (Ile Ife)"))

ggplot(AllelesInt, aes(x = Allele_Count)) +
  geom_histogram(binwidth = 1, color = "black", fill = "darkorchid", alpha = 0.7) +
  labs(
    title = "Distribution of Allele Counts in Villages of Interest",
    x = "Number of Alleles",
    y = "Frequency"
  ) +
  theme_minimal()

# Step 3: Compute summary statistics
summary_int <- AllelesInt %>%
  summarise(
    Median = median(Allele_Count),
    Mean = mean(Allele_Count),
    Min = min(Allele_Count),
    Max = max(Allele_Count),
    Q1 = quantile(Allele_Count, 0.25),
    Q3 = quantile(Allele_Count, 0.75)
  )

# Print the results
summary_int

