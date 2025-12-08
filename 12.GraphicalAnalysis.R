
setwd("C:/Users/Jan_S/Desktop/2024_11_12_Post-Acacia-Final")
setwd("C:/Users/Jan_S/Desktop/2025-01-Mastomys-Final-Really")

{r load-packages, message=FALSE, results='hide'}

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


sum(grepl("^ManaMHC", names(mana_mhc)))


### ALLELES  histograms #########

# Create a new dataframe with just the ManaMHCI columns
mhci_cols <- mana_mhc[, grep("^ManaMHC_", colnames(mana_mhc))]

# Count the number of 1s for each allele
allele_frequencies <- colSums(mhci_cols == 1)


# Convert the allele frequencies into a dataframe for plotting
allele_freq_df <- data.frame(Allele = names(allele_frequencies), Frequency = allele_frequencies)




# Calculate 10% of the total number of samples (527)
total_samples <- nrow(mana_mhc)
threshold <- 0.1 * total_samples




# Create the bar plot with a horizontal line at 10%
ggplot(allele_freq_df, aes(x = Allele, y = Frequency)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +  # Add horizontal line
  labs(x = "Allele", y = "Frequency of Presence (1's)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = nrow(allele_freq_df), y = threshold + 5, label = "10% threshold", color = "red", angle = 0, hjust = 1)




# Create a new dataframe with just the ManaMHCI columns
mhci_cols <- mana_mhc[, grep("^ManaMHC_", colnames(mana_mhc))]

# Count the number of 1s for each allele
allele_frequencies <- colSums(mhci_cols == 1)

# Convert the allele frequencies into a dataframe for plotting
allele_freq_df <- data.frame(Allele = names(allele_frequencies), Frequency = allele_frequencies)

# Reorder Allele levels based on Frequency
allele_freq_df$Allele <- factor(allele_freq_df$Allele, levels = allele_freq_df$Allele[order(allele_freq_df$Frequency, decreasing = TRUE)])

# Calculate 10% of the total number of samples (527)
total_samples <- nrow(mana_mhc)
threshold <- 0.1 * total_samples

# Create the bar plot with a horizontal line at 10%, sorted by frequency
ggplot(allele_freq_df, aes(x = Allele, y = Frequency)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +  # Add horizontal line
  labs(x = "Allele", y = "Frequency of Presence (1's)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = nrow(allele_freq_df), y = threshold + 5, label = "10% threshold", color = "red", angle = 0, hjust = 1)



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
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +  # Add horizontal line
  labs(x = "Allele", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate("text", x = length(unique(allele_summary$Allele)), y = threshold + 5, 
           label = "10% threshold", color = "red", angle = 0, hjust = 1)


################## Per Country stacked

library(dplyr)
library(reshape2)
library(ggplot2)
library(forcats)

# Extract relevant columns (Country, Village, and ManaMHC_ columns)
mhci_cols <- mana_mhc[, c("Country", "Village", grep("^ManaMHC_", colnames(mana_mhc), value = TRUE))]

# Reshape to long format
allele_long <- melt(mhci_cols, id.vars = c("Country", "Village"), 
                    variable.name = "Allele", value.name = "Presence")

# Filter to include only alleles with Presence == 1
allele_long <- allele_long[allele_long$Presence == 1, ]

# Summarize frequencies for each allele by Village within each Country
allele_summary <- allele_long %>%
  group_by(Country, Village, Allele) %>%
  summarise(Frequency = n(), .groups = "drop")

# Create separate plots for each country
plots <- list()
for (country in unique(allele_summary$Country)) {
  # Filter data for the current country
  country_data <- allele_summary %>%
    filter(Country == country)
  
  # Calculate total frequency for alleles in the current country
  allele_order <- country_data %>%
    group_by(Allele) %>%
    summarise(TotalFrequency = sum(Frequency), .groups = "drop") %>%
    arrange(desc(TotalFrequency))
  
  # Reorder alleles for the current country
  country_data <- country_data %>%
    mutate(Allele = factor(Allele, levels = allele_order$Allele))
  
  # Create the plot for the current country
  plot <- ggplot(country_data, aes(x = Allele, y = Frequency, fill = Village)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Allele Frequencies in", country),
         x = "Allele", y = "Frequency", fill = "Village") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Store the plot in the list
  plots[[country]] <- plot
}

# Display the plots

# Display the plots
plots[["Guinea"]]  # Plot for Guinea
plots[["Nigeria"]]  # Plot for Nigeria


################## facet wrap villages


library(dplyr)
library(reshape2)
library(ggplot2)
library(forcats)

# Extract relevant columns (Country, Village, and ManaMHC_ columns)
mhci_cols <- mana_mhc[, c("Country", "Village", grep("^ManaMHC_", colnames(mana_mhc), value = TRUE))]

# Reshape to long format
allele_long <- melt(mhci_cols, id.vars = c("Country", "Village"), 
                    variable.name = "Allele", value.name = "Presence")

# Filter to include only alleles with Presence == 1
allele_long <- allele_long[allele_long$Presence == 1, ]

# Summarize frequencies for each allele by Village within each Country
allele_summary <- allele_long %>%
  group_by(Country, Village, Allele) %>%
  summarise(Frequency = n(), .groups = "drop")

# Function to create a plot for a specific country
create_country_plot <- function(country) {
  # Filter data for the given country
  country_data <- allele_summary %>%
    filter(Country == country)
  
  # Calculate total frequency for alleles in the current country
  allele_order <- country_data %>%
    group_by(Allele) %>%
    summarise(TotalFrequency = sum(Frequency), .groups = "drop") %>%
    arrange(desc(TotalFrequency))
  
  # Reorder alleles for the current country
  country_data <- country_data %>%
    mutate(Allele = factor(Allele, levels = allele_order$Allele))
  
  # Create the plot
  ggplot(country_data, aes(x = Allele, y = Frequency, fill = Village)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Village, scales = "free_y") +  # Facet by Village
    labs(title = paste("Allele Frequencies in", country),
         x = "Allele", y = "Frequency", fill = "Village") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Create plots for Nigeria and Guinea
plot_nigeria <- create_country_plot("Nigeria")
plot_guinea <- create_country_plot("Guinea")

# Display the plots
plot_nigeria
plot_guinea


########## SUPERTYPE HISTO ##################




# Create a new dataframe with just the ManaMHCI columns
ST_cols <- mana_mhc[, grep("^Supertype_", colnames(mana_mhc))]

# Count the number of 1s for each allele
ST_frequencies <- colSums(ST_cols == 1)


# Convert the allele frequencies into a dataframe for plotting
ST_freq_df <- data.frame(ST = names(ST_frequencies), Frequency = ST_frequencies)




# Calculate 10% of the total number of samples (527)
total_samples <- nrow(mana_mhc)
threshold <- 0.1 * total_samples


# Create a new dataframe with just the Supertype columns and Country
ST_cols <- mana_mhc[, c("Country", grep("^Supertype_", colnames(mana_mhc), value = TRUE))]

# Melt the dataframe to a long format for ggplot

ST_long <- melt(ST_cols, id.vars = "Country", variable.name = "Supertype", value.name = "Presence")

# Summarize the data: calculate frequencies for each supertype and country
ST_summary <- ST_long %>%
  dplyr::group_by(Supertype, Country) %>%
  dplyr::summarize(Presence = sum(Presence, na.rm = TRUE), .groups = 'drop')


# Sort the Supertypes using dplyr
ST_summary <- ST_summary %>%
  mutate(
    Supertype = factor(
      Supertype,
      levels = unique(Supertype)[order(as.numeric(gsub("Supertype_", "", unique(Supertype))))]
    )
  )

# Check the sorted order
unique(ST_summary$Supertype)

# Create the bar plot with sorted Supertypes
ggplot(ST_summary, aes(x = Supertype, y = Presence, fill = Country)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Nigeria"="cornflowerblue", "Guinea"="darksalmon"))+
  geom_hline(yintercept = (0.05 * total_samples), linetype = "dashed", color = "black") +  # Add horizontal line
  labs(x = "Supertype", y = "Frequency of Presence", fill = "Country") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  annotate(
    "text",
    x = length(unique(ST_summary$Supertype)),
    y = threshold + 5,
    label = "5% threshold",
    color = "black",
    angle = 0,
    hjust = .9
  )
  



##### SUPERTYPES COUNTRY #######


# Create a new dataframe with just the ManaMHCI columns and the Country column
ST_cols_country <- mana_mhc[, c("Country", grep("^Supertype_", colnames(mana_mhc), value = TRUE))]

# Convert to long format for easier grouping
ST_long <- melt(ST_cols_country, id.vars = "Country", variable.name = "Supertype", value.name = "Presence")

# Summarize the frequency of each supertype for each country
ST_summary <- aggregate(Presence ~ Country + Supertype, data = ST_long, sum)

# Calculate 10% of the total samples for each country
country_sample_counts <- table(mana_mhc$Country)  # Count the number of samples per country
ST_summary$Threshold <- 0.1 * country_sample_counts[ST_summary$Country]

# Calculate the percentage of presence within each country
ST_summary <- ST_summary %>%
  group_by(Country) %>%
  mutate(Percentage = Presence / sum(Presence) * 100) %>%
  ungroup()

# Sort the Supertypes using dplyr
ST_summary <- ST_summary %>%
  mutate(
    Supertype = factor(
      Supertype,
      levels = unique(Supertype)[order(as.numeric(gsub("Supertype_", "", unique(Supertype))))]
    )
  )






# Updated ggplot with manual color scale


ggplot(ST_summary, aes(x = Supertype, y = Percentage, fill = Country)) +
  geom_bar(stat = "identity", position = position_dodge(width = .8), width = .8) +
  labs(x = "Supertype", y = "Percentage of Presence", title = "Supertype Percentages by Country") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)), # Add space below x-axis title
    panel.grid.major.x = element_blank() # Optional: remove vertical grid lines
  ) +
  scale_x_discrete(expand = expansion(add = c(.5, .5))) + # Add space between categories
  scale_fill_manual(values = c("Nigeria" = "cornflowerblue", "Guinea" = "darksalmon", "Other" = "gold"))

################# SUPERTYPES STACKED PER VILLAGE 


# Extract relevant columns (Country and Supertype_ columns)
ST_cols_country <- mana_mhc[, c("Country", "Village", grep("^Supertype_", colnames(mana_mhc), value = TRUE))]

# Reshape to long format
supertype_long <- melt(ST_cols_country, id.vars = c("Country", "Village"), 
                       variable.name = "Supertype", value.name = "Presence")

# Filter to include only supertypes with Presence == 1
supertype_long <- supertype_long[supertype_long$Presence == 1, ]

# Summarize frequencies for each supertype by Village within each Country
supertype_summary <- supertype_long %>%
  group_by(Country, Village, Supertype) %>%
  summarise(Frequency = n(), .groups = "drop")

# Function to create a plot for a specific country
create_supertype_plot <- function(country) {
  # Filter data for the given country
  country_data <- supertype_summary %>%
    filter(Country == country)
  
  # Calculate total frequency for supertypes in the current country
  supertype_order <- country_data %>%
    group_by(Supertype) %>%
    summarise(TotalFrequency = sum(Frequency), .groups = "drop") %>%
    arrange(desc(TotalFrequency))
  
  # Reorder supertypes for the current country
  country_data <- country_data %>%
    mutate(Supertype = factor(Supertype, levels = supertype_order$Supertype))
  
  # Create the stacked bar plot
  ggplot(country_data, aes(x = Supertype, y = Frequency, fill = Village)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste("Supertype Frequencies in", country),
         x = "Supertype", y = "Frequency", fill = "Village") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Create plots for Nigeria and Guinea
plot_nigeria_st <- create_supertype_plot("Nigeria")
plot_guinea_st <- create_supertype_plot("Guinea")

# Display the plots
plot_nigeria_st
plot_guinea_st


################## SUPERTYPES FACET WRAP VILLAGE

# Extract relevant columns (Country, Village, and Supertype_ columns)
ST_cols_country <- mana_mhc[, c("Country", "Village", grep("^Supertype_", colnames(mana_mhc), value = TRUE))]

# Reshape to long format
supertype_long <- melt(ST_cols_country, id.vars = c("Country", "Village"), 
                       variable.name = "Supertype", value.name = "Presence")

# Filter to include only supertypes with Presence == 1
supertype_long <- supertype_long[supertype_long$Presence == 1, ]

# Summarize frequencies for each supertype by Village within each Country
supertype_summary <- supertype_long %>%
  group_by(Country, Village, Supertype) %>%
  summarise(Frequency = n(), .groups = "drop")

# Function to create a facet-wrapped plot for a specific country
create_facet_plot <- function(country) {
  # Filter data for the given country
  country_data <- supertype_summary %>%
    filter(Country == country)
  
  # Calculate total frequency for supertypes in the current country
  supertype_order <- country_data %>%
    group_by(Supertype) %>%
    summarise(TotalFrequency = sum(Frequency), .groups = "drop") %>%
    arrange(desc(TotalFrequency))
  
  # Reorder supertypes for the current country
  country_data <- country_data %>%
    mutate(Supertype = factor(Supertype, levels = supertype_order$Supertype))
  
  # Create the facet-wrapped plot
  ggplot(country_data, aes(x = Supertype, y = Frequency, fill = Village)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Village, scales = "free_y") +  # Facet by Village
    labs(title = paste("Supertype Frequencies in", country, "by Village"),
         x = "Supertype", y = "Frequency", fill = "Village") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Create facet-wrapped plots for Nigeria and Guinea
facet_plot_nigeria <- create_facet_plot("Nigeria")
facet_plot_guinea <- create_facet_plot("Guinea")

# Display the plots
facet_plot_nigeria
facet_plot_guinea


###################### Facet Wrap, fixed y-scale 


# Extract relevant columns (Country, Village, and Supertype_ columns)
ST_cols_country <- mana_mhc[, c("Country", "Village", grep("^Supertype_", colnames(mana_mhc), value = TRUE))]

# Ensure Supertypes are ordered numerically
ST_summary <- ST_cols_country %>%
  mutate(across(starts_with("Supertype_"), as.numeric))  # Ensure numerical transformation if needed

# Reshape to long format
supertype_long <- melt(ST_cols_country, id.vars = c("Country", "Village"), 
                       variable.name = "Supertype", value.name = "Presence")

# Ensure correct supertype ordering by converting to a factor
supertype_long <- supertype_long %>%
  mutate(Supertype = factor(Supertype, levels = paste0("Supertype_", sprintf("%02d", 1:19))))

# Filter to include only supertypes with Presence == 1
supertype_long <- supertype_long[supertype_long$Presence == 1, ]

# Summarize frequencies for each supertype by Village within each Country
supertype_summary <- supertype_long %>%
  group_by(Country, Village, Supertype) %>%
  summarise(Frequency = n(), .groups = "drop")

# Function to create a facet-wrapped plot for a specific country with a fixed y-scale
create_facet_plot_fixed_y <- function(country) {
  country_data <- supertype_summary %>%
    filter(Country == country) %>%
    mutate(Supertype = factor(Supertype, levels = paste0("Supertype_", sprintf("%02d", 1:19))))
  
  # Create the facet-wrapped plot with fixed y-scale
  ggplot(country_data, aes(x = Supertype, y = Frequency, fill = Village)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ Village) +
    ylim(0, 110) +
    labs(title = paste("Supertype Frequencies in", country, "by Village"),
         x = "Supertype", y = "Frequency", fill = "Village") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Create plots for Nigeria and Guinea
facet_plot_nigeria_fixed_y <- create_facet_plot_fixed_y("Nigeria")
facet_plot_guinea_fixed_y <- create_facet_plot_fixed_y("Guinea")

# Display plots
facet_plot_nigeria_fixed_y
facet_plot_guinea_fixed_y

