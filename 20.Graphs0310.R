setwd("C:/Users/Jan_S/Desktop/2025-01-Mastomys-Final-Really")



rm(list=ls())

x <-
  c(
    "plyr", "tidyverse", "vegan", "RColorBrewer", 
    "rstatix", "car", "ggplot2", "vegan", "ggrepel", 
    "DT", "ggpubr", "reshape2", "dplyr", "remotes", 
    "effects", "glmmTMB", "ggeffects", "lme4", "cooccur", "sjPlot")


lapply(x, function(y) {
  # check if installed, if not install
  if (!y %in% installed.packages()[, "Package"])
    install.packages(y)
  
  # load package
  try(require(y, character.only = T), silent = T)
})


#check loaded Rversion for reporting in the paper

# R.Version()
# Plot the frequencies of the alleles using ggplot2
library(ggplot2)

detach("package:dplyr", unload = TRUE) 
library(dplyr)





mana_mhc <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv")


str(mana_mhc) #### check parameters 
### sex, ID, village etc are all characters (could be turned to factors possibly)
#e.g.
mana_mhc$Village<-as.factor(mana_mhc$Village)
mana_mhc$ID<-as.factor(mana_mhc$ID)
mana_mhc$Country<-as.factor(mana_mhc$Country)
mana_mhc$LASV_positive<-as.factor(mana_mhc$ LASV_positive)
mana_mhc$IgG_positive<-as.factor(mana_mhc$IgG_positive)
mana_mhc$State<-as.factor(mana_mhc$State)

# Rename the supertype columns in mana_mhc
colnames(mana_mhc) <- gsub("Supertype_(\\d)$", "Supertype_0\\1", colnames(mana_mhc))

# Verify the change
colnames(mana_mhc)


#mana_mhc <- subset(mana_mhc, mana_mhc$Village %in% c("Damania", "Sokourala", "Ekpoma", "Ebudin", "Owo"))
#allele_columns <- grep("^ManaMHC_", colnames(mana_mhc), value = TRUE)
#length(allele_columns)
#############




# Extract ManaMHCI columns and Country
mhci_cols <- mana_mhc[, c("Country", grep("^ManaMHC_", colnames(mana_mhc), value = TRUE))]



# Reshape data to long format
allele_long <- melt(mhci_cols, id.vars = "Country", variable.name = "Allele", value.name = "Presence")
allele_long <- allele_long %>% filter(!is.na(Presence))



# Summarize frequencies for each allele and country
allele_summary <- allele_long %>%
  filter(Presence == 1) %>%  # Only consider rows where Presence is 1
  group_by(Allele, Country) %>%
  summarise(Frequency = n(), .groups = 'drop')



# Calculate the total frequency of presence for each allele across all countries
allele_total_freq <- allele_summary %>%
  group_by(Allele) %>%
  summarise(TotalFrequency = sum(Frequency), .groups = 'drop') %>%
  arrange(desc(TotalFrequency))  # Sort alleles by total frequency

# Reorder Allele levels based on TotalFrequency
allele_summary$Allele <- factor(allele_summary$Allele, levels = allele_total_freq$Allele)

# Calculate 10% threshold
total_samples <- nrow(mana_mhc)
threshold <- 0.1 * total_samples

nigeria_color <- "#008753"
guinea_color <- "#D19D00"

# Create the stacked bar plot with color coding by country
ggplot(allele_summary, aes(x = Allele, y = Frequency, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria" = nigeria_color, "Guinea" = guinea_color)) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +  # Add horizontal line
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = length(unique(allele_summary$Allele)), y = threshold + 13, 
           label = "10% threshold", color = "black", angle = 0, hjust = 1, size = 3)



# Nigeria cornflowerblue #6495ed
# Guinea #fea600
############# ABOVE/Equal 10 % ###########################

# Identify alleles above the threshold
above10_alleles_list <- allele_total_freq$Allele[allele_total_freq$TotalFrequency >= 66.5] # 66.5 is the value of 10% threshold

# Number of alleles with total frequency >= 10%
length(above10_alleles_list)

# Filter allele_summary based on the selected alleles #     here each allele is listed, seperatly for each country if in both
above10_alleles <- subset(allele_summary, Allele %in% above10_alleles_list)



# Reorder Allele levels correctly
above10_alleles$Allele <- factor(above10_alleles$Allele, levels = above10_alleles_list)


# Create the stacked bar plot with color coding by country
ggplot(above10_alleles, aes(x = Allele, y = Frequency, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria" = nigeria_color, "Guinea" = guinea_color)) +
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) 



############################ Flip Coords


allele_order <- allele_total_freq %>% 
  filter(Allele %in% above10_alleles_list) %>% 
  arrange(TotalFrequency) %>%  # Ascending order first...
  pull(Allele)

# Apply the new order to the factor levels
above10_alleles$Allele <- factor(above10_alleles$Allele, levels = allele_order)

# Create the stacked bar plot with color coding by country and flipped axes
ggplot(above10_alleles, aes(x = Allele, y = Frequency, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria" = nigeria_color, "Guinea" = guinea_color)) +
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) +  # Adjust for readability
  coord_flip()  # Flip axes so highest frequency is on top


###################             with *         ######################
# Create new versions of the dataframes with renamed alleles
allele_total_freq2 <- allele_total_freq
above10_alleles2 <- above10_alleles

# Replace "ManaMHC_" with "*" in both dataframes
allele_total_freq2$Allele <- gsub("ManaMHC_", "*", allele_total_freq2$Allele)
above10_alleles2$Allele <- gsub("ManaMHC_", "*", above10_alleles2$Allele)

# ALSO update above10_alleles_list so filtering works!
above10_alleles_list2 <- gsub("ManaMHC_", "*", above10_alleles_list)

# Create the new allele order
allele_order2 <- allele_total_freq2 %>% 
  filter(Allele %in% above10_alleles_list2) %>% 
  arrange(TotalFrequency) %>%  # Ascending order first
  pull(Allele)

# Apply the new order to the factor levels
above10_alleles2$Allele <- factor(above10_alleles2$Allele, levels = allele_order2)

# Create the stacked bar plot with color coding by country and flipped axes
ggplot(above10_alleles2, aes(x = Allele, y = Frequency, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria" = nigeria_color, "Guinea" = guinea_color)) +
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10)) +  # Adjust for readability
  coord_flip()  # Flip axes so highest frequency is on top


################## Frequency above/equal 100  (99,75 = 15%)###########################

# Identify alleles above the threshold
above100_alleles_list <- allele_total_freq$Allele[allele_total_freq$TotalFrequency >= 100] 

# Filter allele_summary based on the selected alleles
above100_alleles <- subset(allele_summary, Allele %in% above100_alleles_list)

##### 86 Alleles


# Reorder Allele levels correctly
above100_alleles$Allele <- factor(above100_alleles$Allele, levels = above100_alleles_list)





# Create the stacked bar plot with color coding by country
ggplot(above100_alleles, aes(x = Allele, y = Frequency, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria" = nigeria_color, "Guinea" = guinea_color)) +
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



################################################


# Filter mana_mhc for specific villages
filtered_mana_mhc <- subset(mana_mhc, Village %in% c("Damania", "Sokourala", "Ekpoma", "Ebudin", "Owo"))

# Extract ManaMHC columns and Country
mhci_cols <- mana_mhc[, c("Country", "Village", grep("^ManaMHC_", colnames(mana_mhc), value = TRUE))]

# Reshape to long format
allele_long <- melt(mhci_cols, id.vars = c("Country", "Village"), variable.name = "Allele", value.name = "Presence")

allele_long <- allele_long %>% filter(!is.na(Presence))

# Step 4: Summarize total allele frequencies (global, across all samples)
allele_summary <- allele_long %>%
  filter(Presence == 1) %>%
  group_by(Allele) %>%
  summarise(TotalFrequency = n(), .groups = 'drop') %>%
  arrange(desc(TotalFrequency))

# Step 5: Apply 10% threshold (based on full sample size)
threshold <- 0.1 * nrow(mana_mhc)
above10_alleles_list <- allele_summary %>%
  filter(TotalFrequency >= threshold) %>%
  pull(Allele)

# Apply threshold across all alleles

allele_long_filtered <- allele_long %>%
  filter(Allele %in% above10_alleles_list)

target_villages <- c("Damania", "Sokourala", "Ekpoma", "Ebudin", "Owo")
allele_long_villages <- allele_long_filtered %>%
  filter(Village %in% target_villages)

# Step 8 (Optional): Check how many alleles are retained in the villages
alleles_in_villages <- unique(allele_long_villages$Allele)
length(alleles_in_villages)


###########




