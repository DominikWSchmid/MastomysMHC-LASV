setwd("C:/Users/jansa/Desktop/AYO/2025-Mastomys")

load("19.1PERMANOVA_Guinea.RData")
W <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv")


str(W)

# Load necessary package
library(vegan)
library(ggplot2)

# Filter data for Guinea
Guinea_data <- W[W$Country == "Guinea", ]
Guinea_data <- na.omit(Guinea_data)
# Extract allele composition matrix (response variable)
gallele_matrix <- Guinea_data[, grep("^ManaMHC_", names(Guinea_data))]

# Ensure response variable is numeric
gallele_matrix <- as.matrix(gallele_matrix)

# Create explanatory variable (Village)
gvillage_factor <- factor(Guinea_data$Village)

############################ PERMANOVA JACCARD 
# Compute Jaccard distance matrix
gjaccard_dist <- vegdist(gallele_matrix, method = "jaccard", binary = TRUE)

# Run PERMANOVA with Jaccard distance
set.seed(123) # For reproducibility
gpermanova_jaccard <- adonis2(gjaccard_dist ~ gvillage_factor, data = Guinea_data)

# View results
print(gpermanova_jaccard)

#############################  PERMANOVA EUCLIDEAN

# Compute Euclidean distance matrix
geuclidean_dist <- dist(gallele_matrix, method = "euclidean")

# Run PERMANOVA with Euclidean distance
set.seed(123) # For reproducibility
gpermanova_euclidean <- adonis2(geuclidean_dist ~ gvillage_factor, data = Guinea_data)

# View results
print(gpermanova_euclidean)


#############################  PERMANOVA Manhattan

# Compute manhattan distance matrix
gmanhattan_dist <- dist(gallele_matrix, method = "manhattan")

# Run PERMANOVA with manhattan distance
set.seed(123) # For reproducibility
gpermanova_manhattan <- adonis2(gmanhattan_dist ~ gvillage_factor, data = Guinea_data)

# View results
print(gpermanova_manhattan)


################################


# Compute Manhattan distance matrix
gman_dist <- vegdist(gallele_matrix, method = "manhattan")

# Perform PCoA
gpcoa_result <- cmdscale(gman_dist, eig = TRUE, k = 2) # k = 2 for 2D visualization

# Extract coordinates for plotting
gpcoa_coords <- as.data.frame(gpcoa_result$points)
colnames(gpcoa_coords) <- c("PC1", "PC2")

# Add village information for visualization
gpcoa_coords$Village <- Guinea_data$Village

# Plot PCoA results
ggplot(gpcoa_coords, aes(x = PC1, y = PC2, color = Village)) +
  stat_ellipse(aes(fill = Village), type = "norm", level = 0.95, alpha = 0.2, geom = "polygon") + # 95% confidence ellipsoids
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCoA of Allele Composition (Manhattan Distance)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right")



########### PAIRWISE ############

library(pairwiseAdonis)

# Jaccard Distance Pairwise Comparisons
pairwise_jaccard <- pairwise.adonis(gjaccard_dist, factors = gvillage_factor, sim.function = vegdist, sim.method = "jaccard", p.adjust.m = "bonferroni")
print(pairwise_jaccard)

# Euclidean Distance Pairwise Comparisons
pairwise_euclidean <- pairwise.adonis(geuclidean_dist, factors = gvillage_factor, sim.function = dist, sim.method = "euclidean", p.adjust.m = "bonferroni")
print(pairwise_euclidean)

# Manhattan Distance Pairwise Comparisons
pairwise_manhattan <- pairwise.adonis(gmanhattan_dist, factors = gvillage_factor, sim.function = dist, sim.method = "manhattan", p.adjust.m = "bonferroni")
print(pairwise_manhattan)
