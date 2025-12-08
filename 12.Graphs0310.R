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


# Create the stacked bar plot with color coding by country
ggplot(allele_summary, aes(x = Allele, y = Frequency, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria" = "cornflowerblue", "Guinea" = "darksalmon")) +
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

# Filter allele_summary based on the selected alleles
above10_alleles <- subset(allele_summary, Allele %in% above10_alleles_list)


#### 138 Alleles


# Reorder Allele levels correctly
above10_alleles$Allele <- factor(above10_alleles$Allele, levels = above10_alleles_list)


# Create the stacked bar plot with color coding by country
ggplot(above10_alleles, aes(x = Allele, y = Frequency, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria" = "cornflowerblue", "Guinea" = "darksalmon")) +
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) 

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
  scale_fill_manual(values = c("Nigeria" = "cornflowerblue", "Guinea" = "darksalmon")) +
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



################################################








