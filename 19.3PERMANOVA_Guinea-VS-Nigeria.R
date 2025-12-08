
# Load necessary package
library(vegan)
library(ggplot2)



load("19.3PERMANOVA_Guina-VS-Nigeria.RData")

W <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv")

str(W)
africa_data <- W

africa_data <- africa_data[africa_data$Village %in% c("Ebudin", "Ekpoma", "Owo","Aba gboro (Ile Ife)", "Damania", "Sokourala", "Tambaya", "Sonkonia"), ]
allele_matrix <- africa_data[, grep("^ManaMHC_", names(africa_data))]
supertype_matrix <- africa_data[, grep("^Supertype_", names(africa_data))]


allele_matrix <- as.matrix(allele_matrix)

# Create explanatory variable (country)
country_factor <- factor(africa_data$Country)
village_factor <- factor(africa_data$Village)

############################ PERMANOVA JACCARD 
# Compute Jaccard distance matrix
jaccard_dist <- vegdist(allele_matrix, method = "jaccard", binary = TRUE)

# Run PERMANOVA with Jaccard distance
set.seed(123) # For reproducibility
permanova_jaccard <- adonis2(jaccard_dist ~ country_factor, data = africa_data)

# View results
print(permanova_jaccard)

#############################  PERMANOVA EUCLIDEAN

# Compute Euclidean distance matrix
euclidean_dist <- dist(allele_matrix, method = "euclidean")

# Run PERMANOVA with Euclidean distance
set.seed(123) # For reproducibility
permanova_euclidean <- adonis2(euclidean_dist ~ country_factor, data = africa_data)

# View results
print(permanova_euclidean)


#############################  PERMANOVA Manhattan

# Compute manhattan distance matrix
manhattan_dist <- dist(allele_matrix, method = "manhattan")
manhattan_dist_sup <- dist(supertype_matrix, method = "manhattan")

# Run PERMANOVA with manhattan distance
set.seed(123) # For reproducibility
permanova_manhattan <- adonis2(manhattan_dist ~ country_factor, data = africa_data)
permanova_manhattan_sup <- adonis2(manhattan_dist_sup ~ country_factor, data = africa_data)
permanova_manhattan <- adonis2(manhattan_dist ~ village_factor, data = africa_data)
permanova_manhattan_sup <- adonis2(manhattan_dist_sup ~ village_factor, data = africa_data)

# View results
print(permanova_manhattan)
print(permanova_manhattan_sup)

print(permanova_euclidean)

print(permanova_jaccard)
################################


# Compute Manhattan distance matrix
man_dist <- vegdist(allele_matrix, method = "manhattan")

# Perform PCoA
pcoa_result <- cmdscale(man_dist, eig = TRUE, k = 2) # k = 2 for 2D visualization

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("PC1", "PC2")

# Add country information for visualization
pcoa_coords$Country <- africa_data$Country

# Base colors for Nigeria and Guinea
nigeria_color <- "#008753"
guinea_color <- "#D19D00"
country_colors <- c("Nigeria" = nigeria_color, "Guinea" = guinea_color)


# Plot PCoA results
ggplot(pcoa_coords, aes(x = PC1, y = PC2, color = Country)) +
  stat_ellipse(aes(fill = Country), type = "norm", level = 0.95, alpha = 0.2, geom = "polygon") + # 95% confidence ellipsoids
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCoA of Allele Composition (Manhattan Distance)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right")+
  scale_color_manual(values = country_colors) + # Apply the country colors
  scale_fill_manual(values = country_colors)   # Apply the country colors to the fill









######################PAIRWISE Comparison############################
library(pairwiseAdonis)


# Jaccard Distance Pairwise Comparisons
pairwise_jaccard <- pairwise.adonis(jaccard_dist, factors = country_factor, sim.function = vegdist, sim.method = "jaccard", p.adjust.m = "bonferroni")
print(pairwise_jaccard)

# Euclidean Distance Pairwise Comparisons
pairwise_euclidean <- pairwise.adonis( euclidean_dist, factors =  country_factor, sim.function = dist, sim.method = "euclidean", p.adjust.m = "bonferroni")
print(pairwise_euclidean)

# Manhattan Distance Pairwise Comparisons
pairwise_manhattan <- pairwise.adonis(manhattan_dist, factors = village_factor, sim.function = dist, sim.method = "manhattan", p.adjust.m = "bonferroni")
print(pairwise_manhattan)


######################### ANOSIM   ###############


# Run ANOSIM with different distance metrics
set.seed(123)  # For reproducibility
anosim_manhattan <- anosim(vegdist(allele_matrix, method = "manhattan"), africa_data$Village)
anosim_euclidean <- anosim(vegdist(allele_matrix, method = "euclidean"), africa_data$Village)
anosim_jaccard <- anosim(vegdist(allele_matrix, method = "jaccard"), africa_data$Village)

# Print results
print(anosim_manhattan)
print(anosim_euclidean)
print(anosim_jaccard)


################ NMDS #############

library(ggrepel)
library(dplyr)

# Compute NMDS for different distance metrics
#nmds_manhattan <- metaMDS(allele_matrix, distance = "manhattan", k = 3, trymax = 100)
#nmds_euclidean <- metaMDS(allele_matrix, distance = "euclidean", k = 3, trymax = 100)
#nmds_jaccard <- metaMDS(allele_matrix, distance = "jaccard", k = 3, trymax = 100)




#saveRDS(nmds_manhattan, file = "nmds_manhattan_k=3.rds")
#saveRDS(nmds_euclidean, file = "nmds_euclidean_k=3.rds")
#saveRDS(nmds_jaccard, file = "nmds_jaccard_k=3.rds")


nmds_manhattan <- readRDS("nmds_manhattan_k=3.rds")
nmds_euclidean <- readRDS("nmds_euclidean_k=3.rds")
nmds_jaccard <- readRDS("nmds_jaccard_k=3.rds")


# Extract NMDS scores and add village info
nmds_data_ja <- data.frame(
  NMDS1 = nmds_jaccard$points[, 1],
  NMDS2 = nmds_jaccard$points[, 2],
  Village = africa_data$Village,
  Country = africa_data$Country
)


# Extract NMDS scores and add village info
nmds_data_ma <- data.frame(
  NMDS1 = nmds_manhattan$points[, 1],
  NMDS2 = nmds_manhattan$points[, 2],
  Village = africa_data$Village,
  Country = africa_data$Country
)


# Extract NMDS scores and add village info
nmds_data_eu <- data.frame(
  NMDS1 = nmds_euclidean$points[, 1],
  NMDS2 = nmds_euclidean$points[, 2],
  Village = africa_data$Village,
  Country = africa_data$Country
)


# Define colors for villages
village_colors <- c(
  "Aba gboro (Ile Ife)" = "#0EFDA2",
  "Ebudin" = "#03C97D",
  "Ekpoma" = "#038A56",
  "Owo" = "#024B2F",
  "Damania" = "#4F3C02",
  "Sokourala" = "#8F6C03",
  "Sonkonia" = "#CE9C03",
  "Tambaya" = "#FEC313"
)

# Plot NMDS-jaccard
nmds_plot_ja <- ggplot(nmds_data_ja, aes(x = NMDS1, y = NMDS2, color = Village)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
  # geom_text_repel(aes(label = Village), size = 3, max.overlaps = 100) +
  theme_minimal() +
  scale_color_manual(values = village_colors) +
  ggtitle("NMDS of Allele Composition (Jaccard Distance)") +
  theme(legend.position = "right")


# Plot NMDS-manhattan
nmds_plot_ma <- ggplot(nmds_data_ma, aes(x = NMDS1, y = NMDS2, color = Village)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = 5) +
#  geom_text_repel(aes(label = Village), size = 3, max.overlaps = 100) +
  theme_minimal() +
  scale_color_manual(values = village_colors) +
  ggtitle("NMDS of Allele Composition (Manhattan Distance)") +
  theme(legend.position = "right")

nmds_plot_ma

# Plot NMDS_
nmds_plot_eu <- ggplot(nmds_data_eu, aes(x = NMDS1, y = NMDS2, color = Village)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 2) +
#  geom_text_repel(aes(label = Village), size = 3, max.overlaps = 100) +
  theme_minimal() +
  scale_color_manual(values = village_colors) +
  ggtitle("NMDS of Allele Composition (Euclidean Distance)") +
  theme(legend.position = "right")



# Print NMDS stress value
stress_value_ja <- nmds_jaccard$stress
print(paste("NMDS Stress:", round(stress_value_ja, 3)))

# Print NMDS stress value
stress_value_ma <- nmds_manhattan$stress
print(paste("NMDS Stress:", round(stress_value_ma, 3)))

# Print NMDS stress value
stress_value_eu <- nmds_euclidean$stress
print(paste("NMDS Stress:", round(stress_value_eu, 3)))


# Show plot
print(nmds_plot_ja)
print(nmds_plot_ma)
print(nmds_plot_eu)

###############


# Plot NMDS-manhattan
nmds_plot_ma <- ggplot(nmds_data_ma, aes(x = NMDS1, y = NMDS2, color = Village)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 6) +
  #  geom_text_repel(aes(label = Village), size = 3, max.overlaps = 100) +
  theme_minimal() +
  scale_color_manual(values = village_colors) +
  ggtitle("NMDS of Allele Composition (Manhattan Distance)") +
  theme(legend.position = "right")


print(nmds_plot_ma)

#######



# 1. Create the VillageCountry factor
nmds_data_ma$VillageCountry <- paste(nmds_data_ma$Country, "-", nmds_data_ma$Village)

# 2. Order the levels of the VillageCountry factor
nmds_data_ma <- within(nmds_data_ma,
                       VillageCountry <- factor(VillageCountry,
                                                levels = unique(paste(Country, "-", Village)[order(Country, Village)])))

# 3. Create a new color vector mapping VillageCountry to colors
village_country_levels <- levels(nmds_data_ma$VillageCountry)
village_colors_for_legend <- character(length(village_country_levels))
names(village_colors_for_legend) <- village_country_levels

for (i in seq_along(village_country_levels)) {
  village_name <- gsub("^Nigeria - |^Guinea - ", "", village_country_levels[i])
  if (village_name %in% names(village_colors)) {
    village_colors_for_legend[village_country_levels[i]] <- village_colors[village_name]
  }
}

# 4. Plot NMDS-manhattan with VillageCountry for color and the new color scale
nmds_plot_ma <- ggplot(nmds_data_ma, aes(x = NMDS1, y = NMDS2, color = VillageCountry)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 6) +
  theme_minimal() +
  scale_color_manual(values = village_colors_for_legend,
                     labels = levels(nmds_data_ma$VillageCountry)) +
  ggtitle("NMDS of Allele Composition (Manhattan Distance)") +
  labs(color = "Villages") +
  theme(legend.position = "right")

# Print the plot
print(nmds_plot_ma)
