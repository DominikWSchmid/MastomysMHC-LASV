setwd("C:/Users/Jan_S/Desktop/2024_11_12_Post-Acacia-Final")


mergedG <- read.csv("Nigeria_Guinea_Wide_ST.csv")

######################### Figures #################

library(lubridate)
library(dplyr)
library(ggplot2)
library(tidyr)

# Step 1: Count unique IDs per date
date_counts <- mergedG %>%
  filter(!is.na(Date_Capture)) %>% # Remove rows without dates
  group_by(Date_Capture) %>%      # Group by capture date
  summarise(ID_count = n_distinct(ID)) # Count unique IDs

# Step 2: Plot

ggplot(date_counts, aes(x = as.Date(Date_Capture), y = ID_count)) +
  geom_line(color = "#c37992", size = 1) +   # Line for trends
  geom_point(color = "red", size = 2) +  # Points for individual dates
  labs(
    title = "Number of Unique IDs Over Time",
    x = "Date",
    y = "Number of Unique IDs"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
################## Monthly ####################

# Step 1: Group date points by month
date_counts_monthly <- mergedG %>%
  filter(!is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month) %>%
  summarise(ID_count = n_distinct(ID))

variant_counts_monthly <- long_data %>%
  filter(ManaMHCI %in% c("ManaMHCI.006", "ManaMHCI.008", "ManaMHCI.021")) %>%
  filter(!is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month, ManaMHCI) %>%
  summarise(Variant_count = n_distinct(ID), .groups = "drop")

# Step 2: Plotting
ggplot() +
  # Add smoothed line for all IDs by month
  geom_smooth(data = date_counts_monthly, aes(x = Month, y = ID_count), 
              method = "loess", se = FALSE, color = "grey", size = 3, alpha = 0.7) +
  geom_point(data = date_counts_monthly, aes(x = Month, y = ID_count), 
             color = "black", size = 2) +
  # Overlay smoothed lines for specific ManaMHCI variants by month
  geom_smooth(data = variant_counts_monthly, aes(x = Month, y = Variant_count, 
                                                 color = ManaMHCI, group = ManaMHCI),
              method = "loess", se = FALSE, size = 1.2) +
  geom_point(data = variant_counts_monthly, aes(x = Month, y = Variant_count, 
                                                color = ManaMHCI), size = 2) +
  # Add labels and a legend
  labs(
    title = "Monthly Trends of Unique IDs and ManaMHCI Variants",
    x = "Month",
    y = "Count",
    color = "ManaMHCI Variant"
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + # Adjust X-axis for readability
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#
######## Relative frequency ####################

# Step 1: Group by Month and Calculate Total IDs
total_ids_monthly <- mergedG %>%
  filter(!is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month) %>%
  summarise(Total_IDs = n_distinct(ID), .groups = "drop")

# Step 2: Group by Month and Allele, and Calculate Counts
allele_counts_monthly <- long_data %>%
  filter(ManaMHCI %in% c("ManaMHCI.006", "ManaMHCI.008", "ManaMHCI.021")) %>%
  filter(!is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month, ManaMHCI) %>%
  summarise(Allele_Count = n_distinct(ID), .groups = "drop")

# Step 3: Merge and Calculate Proportions
allele_frequencies <- merge(allele_counts_monthly, total_ids_monthly, by = "Month") %>%
  mutate(Proportion = Allele_Count / Total_IDs)

# Step 4: Plot Proportions Over Time
ggplot() +
  # Smoothed lines for proportions of selected alleles
  geom_smooth(data = allele_frequencies, aes(x = Month, y = Proportion, 
                                             color = ManaMHCI, group = ManaMHCI), 
              method = "loess", se = FALSE, size = 1.2) +
  geom_point(data = allele_frequencies, aes(x = Month, y = Proportion, color = ManaMHCI), 
             size = 2) +
  # Add labels and a legend
  labs(
    title = "Relative Frequency of Selected ManaMHCI Alleles Over Time",
    x = "Month",
    y = "Proportion of IDs",
    color = "ManaMHCI Variant"
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + # Adjust X-axis for readability
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



############ plus positive #################

# Step 1: Calculate LASV-positive Proportions
lasv_positive_monthly <- mergedG %>%
  filter(!is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month) %>%
  summarise(
    LASV_Positive_Count = sum(LASV_positive == "pos", na.rm = TRUE),
    Total_IDs = n_distinct(ID),
    LASV_Proportion = LASV_Positive_Count / Total_IDs,
    .groups = "drop"
  )

# Step 2: Combine LASV Proportions with Allele Data
allele_frequencies_lasv <- merge(
  allele_frequencies, lasv_positive_monthly[, c("Month", "LASV_Proportion")],
  by = "Month", all.x = TRUE
)

# Step 3: Plot Proportions of Alleles and LASV Positivity
ggplot() +
  # Smoothed lines for ManaMHCI allele proportions
  geom_smooth(data = allele_frequencies_lasv, aes(x = Month, y = Proportion, 
                                                  color = ManaMHCI, group = ManaMHCI), 
              method = "loess", se = FALSE, size = 1.2) +
  geom_point(data = allele_frequencies_lasv, aes(x = Month, y = Proportion, color = ManaMHCI), 
             size = 2) +
  # Smoothed line for LASV-positive proportion
  geom_smooth(data = lasv_positive_monthly, aes(x = Month, y = LASV_Proportion), 
              method = "loess", se = FALSE, color = "black", size = 1.2, linetype = "dashed") +
  geom_point(data = lasv_positive_monthly, aes(x = Month, y = LASV_Proportion), 
             color = "black", size = 2) +
  # Add labels and a legend
  labs(
    title = "Relative Frequency of ManaMHCI Alleles and LASV Positivity Over Time",
    x = "Month",
    y = "Proportion of IDs",
    color = "ManaMHCI Variant",
    linetype = ""
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + # Adjust X-axis for readability
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("red", "blue", "green", "black"))




############ only positive villages ################

# Step 1: Filter the data for the specified villages
filtered_data <- mergedG %>%
  filter(Village %in% c("Ebudin", "Ekpoma", "Owo", "Damania", "Sokourala")) %>%
  filter(!is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month"))

# Step 2: Calculate total IDs and LASV-positive proportions per month
total_ids_monthly <- filtered_data %>%
  group_by(Month) %>%
  summarise(
    Total_IDs = n_distinct(ID),
    LASV_Positive_Count = sum(LASV_positive == "pos", na.rm = TRUE),
    LASV_Proportion = LASV_Positive_Count / Total_IDs,
    .groups = "drop"
  )

# Step 3: Calculate allele counts per month
allele_counts_monthly <- long_data %>%
  filter(ManaMHCI %in% c("ManaMHCI.006", "ManaMHCI.008", "ManaMHCI.021")) %>%
  filter(Village %in% c("Ebudin", "Ekpoma", "Owo", "Damania", "Sokourala")) %>%
  filter(!is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month, ManaMHCI) %>%
  summarise(Allele_Count = n_distinct(ID), .groups = "drop")

# Step 4: Merge and calculate allele proportions
allele_frequencies <- merge(allele_counts_monthly, total_ids_monthly, by = "Month") %>%
  mutate(Proportion = Allele_Count / Total_IDs)

# Step 5: Plot the data
ggplot() +
  # Smoothed lines for ManaMHCI allele proportions
  geom_smooth(data = allele_frequencies, aes(x = Month, y = Proportion, 
                                             color = ManaMHCI, group = ManaMHCI), 
              method = "loess", se = FALSE, size = 1.2) +
  geom_point(data = allele_frequencies, aes(x = Month, y = Proportion, color = ManaMHCI), 
             size = 2) +
  # Smoothed line for LASV-positive proportion
  geom_smooth(data = total_ids_monthly, aes(x = Month, y = LASV_Proportion), 
              method = "loess", se = FALSE, color = "black", size = 1.2, linetype = "dashed") +
  geom_point(data = total_ids_monthly, aes(x = Month, y = LASV_Proportion), 
             color = "black", size = 2) +
  # Add labels and a legend
  labs(
    title = "Relative Frequency of ManaMHCI Alleles and LASV Positivity Over Time",
    x = "Month",
    y = "Proportion of IDs",
    color = "ManaMHCI Variant",
    linetype = ""
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") + # Adjust X-axis for readability
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values = c("red", "blue", "green", "black"))



########### two plots ###############
library(gridExtra)

villages_of_interest <- c("Ebudin", "Ekpoma", "Owo", "Damania", "Sokourala")

filtered_data <- filtered_data %>% filter(Village == villages_of_interest)

# Step 1: Filter data for Guinea and Nigeria separately
guinea_data <- filtered_data %>% filter(Country == "Guinea")
nigeria_data <- filtered_data %>% filter(Country == "Nigeria")

# Step 2: Calculate total IDs and LASV-positive proportions for Guinea
guinea_ids_monthly <- guinea_data %>%
  group_by(Month) %>%
  summarise(
    Total_IDs = n_distinct(ID),
    LASV_Positive_Count = sum(LASV_positive == "pos", na.rm = TRUE),
    LASV_Proportion = LASV_Positive_Count / Total_IDs,
    .groups = "drop"
  )

guinea_allele_counts <- long_data %>%
  filter(ManaMHCI %in% c("ManaMHCI.006", "ManaMHCI.008", "ManaMHCI.021")) %>%
  filter(Country == "Guinea", !is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month, ManaMHCI) %>%
  summarise(Allele_Count = n_distinct(ID), .groups = "drop")

guinea_allele_frequencies <- merge(guinea_allele_counts, guinea_ids_monthly, by = "Month") %>%
  mutate(Proportion = Allele_Count / Total_IDs)

# Step 3: Calculate total IDs and LASV-positive proportions for Nigeria
nigeria_ids_monthly <- nigeria_data %>%
  group_by(Month) %>%
  summarise(
    Total_IDs = n_distinct(ID),
    LASV_Positive_Count = sum(LASV_positive == "pos", na.rm = TRUE),
    LASV_Proportion = LASV_Positive_Count / Total_IDs,
    .groups = "drop"
  )

nigeria_allele_counts <- long_data %>%
  filter(ManaMHCI %in% c("ManaMHCI.006", "ManaMHCI.008", "ManaMHCI.021")) %>%
  filter(Country == "Nigeria", !is.na(Date_Capture)) %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month, ManaMHCI) %>%
  summarise(Allele_Count = n_distinct(ID), .groups = "drop")

nigeria_allele_frequencies <- merge(nigeria_allele_counts, nigeria_ids_monthly, by = "Month") %>%
  mutate(Proportion = Allele_Count / Total_IDs)

# Step 4: Create separate plots for Guinea and Nigeria
guinea_plot <- ggplot() +
  
  geom_point(data = guinea_allele_frequencies, aes(x = Month, y = Proportion, color = ManaMHCI), 
             size = 2) +
  geom_line(data = guinea_ids_monthly, aes(x = Month, y = LASV_Proportion), 
           #   method = "loess", se = T, 
           color = "black", size = 1.2, linetype = "dashed") +
  geom_point(data = guinea_ids_monthly, aes(x = Month, y = LASV_Proportion), 
             color = "black", size = 2) +
  geom_line(data = guinea_allele_frequencies, aes(x = Month, y = Proportion, , 
                                                  #     method = "loess", se = F,
                                                  color = ManaMHCI, group = ManaMHCI), size = 1.2) +
  labs(
    title = "Guinea: Relative Frequency of ManaMHCI Alleles and LASV Positivity Over Time",
    x = "Month",
    y = "Proportion of IDs",
    color = "ManaMHCI Variant",
    linetype = ""
  ) +
  theme_minimal()

nigeria_plot <- ggplot() +
  geom_line(data = nigeria_allele_frequencies, aes(x = Month, y = Proportion, 
                                                     color = ManaMHCI, group = ManaMHCI), 
         #     method = "loess", se = F,
            size = 1.2) +
  geom_point(data = nigeria_allele_frequencies, aes(x = Month, y = Proportion, color = ManaMHCI), 
             size = 2) +
  geom_smooth(data = nigeria_ids_monthly, aes(x = Month, y = LASV_Proportion), 
              method = "loess", se = F, color = "black", size = 1.2, linetype = "dashed") +
  geom_point(data = nigeria_ids_monthly, aes(x = Month, y = LASV_Proportion), 
             color = "black", size = 2) +
  labs(
    title = "Nigeria: Relative Frequency of ManaMHCI Alleles and LASV Positivity Over Time",
    x = "Month",
    y = "Proportion of IDs",
    color = "ManaMHCI Variant",
    linetype = ""
  ) +
  theme_minimal()

# Display the plots

grid.arrange(guinea_plot, nigeria_plot, ncol = 1)




####################
# Step 1: Filter the data for Guinea
guinea_data <- mergedG %>%
  filter(Country == "Guinea")

# Step 2: Extract the month and year from Date_Capture and calculate the count of IDs per date and village
guinea_capture_table <- guinea_data %>%
  mutate(Month = floor_date(as.Date(Date_Capture), "month")) %>%
  group_by(Month, Village) %>%
  summarise(ID_Count = n_distinct(ID), .groups = "drop")

# Step 3: Pivot the table to have each village as a separate column
guinea_capture_table_pivot <- guinea_capture_table %>%
  pivot_wider(names_from = Village, values_from = ID_Count, values_fill = list(ID_Count = 0))

# Step 4: Print the table
print(guinea_capture_table_pivot)

#####################

# Define the villages of interest
villages_of_interest <- c("Ebudin", "Ekpoma", "Owo", "Damania", "Sokourala")

# Filter the merged data for Guinea and the selected villages
guinea_data_filtered <- mergedG %>%
  filter(Country == "Guinea", Village %in% villages_of_interest)

# Filter the merged data for Nigeria and the selected villages
nigeria_data_filtered <- mergedG %>%
  filter(Country == "Nigeria", Village %in% villages_of_interest)

# Plot for Guinea
guinea_plot <- ggplot(guinea_data_filtered, aes(x = Date_Capture, color = Village)) +
  geom_line(stat = "count", aes(group = 1), size = 1) +
  labs(title = "ID Capture Over Time in Guinea", x = "Date", y = "Number of IDs") +
  theme_minimal()

# Plot for Nigeria
nigeria_plot <- ggplot(nigeria_data_filtered, aes(x = Date_Capture, color = Village)) +
  geom_line(stat = "count", aes(group = 1), size = 1) +
  labs(title = "ID Capture Over Time in Nigeria", x = "Date", y = "Number of IDs") +
  theme_minimal()

# Display the plots
guinea_plot
nigeria_plot


