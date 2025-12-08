setwd("C:/Users/jansa/Desktop/Ayo/2024_11_12_Post-Acacia-Final")
setwd("C:/Users/Jan_S/Desktop/2025-Mastomys")
setwd("C:/Users/Jan_S/Desktop/2025-01-Mastomys-Final-Really")


load("19.2PERMANOVA_Nigeria.RData")

W <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv")

str(W)

# Load necessary package
library(vegan)
library(ggplot2)

# Filter data for Nigeria
nigeria_data <- W[W$Country == "Nigeria", ]



table(nigeria_data$Village)
nigeria_data <- nigeria_data[nigeria_data$Village %in% c("Ebudin", "Ekpoma", "Owo","Aba gboro (Ile Ife)", "Damania", "Sokourala"), ]

# Extract allele composition matrix (response variable)
allele_matrix <- nigeria_data[, grep("^ManaMHC_", names(nigeria_data))]

# Ensure response variable is numeric
allele_matrix <- as.matrix(allele_matrix)

# Create explanatory variable (Village)
village_factor <- factor(nigeria_data$Village)

############################ PERMANOVA JACCARD 
# Compute Jaccard distance matrix
jaccard_dist <- vegdist(allele_matrix, method = "jaccard", binary = TRUE)

# Run PERMANOVA with Jaccard distance
set.seed(123) # For reproducibility
permanova_jaccard <- adonis2(jaccard_dist ~ village_factor, data = nigeria_data)

# View results
print(permanova_jaccard)

#############################  PERMANOVA EUCLIDEAN

# Compute Euclidean distance matrix
euclidean_dist <- dist(allele_matrix, method = "euclidean")

# Run PERMANOVA with Euclidean distance
set.seed(123) # For reproducibility
permanova_euclidean <- adonis2(euclidean_dist ~ village_factor, data = nigeria_data)

# View results
print(permanova_euclidean)


#############################  PERMANOVA Manhattan

# Compute manhattan distance matrix
manhattan_dist <- dist(allele_matrix, method = "manhattan")

# Run PERMANOVA with manhattan distance
set.seed(123) # For reproducibility
permanova_manhattan <- adonis2(manhattan_dist ~ village_factor, data = nigeria_data)

# View results
print(permanova_manhattan)


################################


# Compute Manhattan distance matrix
man_dist <- vegdist(allele_matrix, method = "manhattan")

# Perform PCoA
pcoa_result <- cmdscale(man_dist, eig = TRUE, k = 2) # k = 2 for 2D visualization

# Extract coordinates for plotting
pcoa_coords <- as.data.frame(pcoa_result$points)
colnames(pcoa_coords) <- c("PC1", "PC2")

# Add village information for visualization
pcoa_coords$Village <- nigeria_data$Village

# Plot PCoA results
ggplot(pcoa_coords, aes(x = PC1, y = PC2, color = Village)) +
  stat_ellipse(aes(fill = Village), type = "norm", level = 0.95, alpha = 0.2, geom = "polygon") + # 95% confidence ellipsoids
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCoA of Allele Composition (Manhattan Distance)",
       x = "PC1",
       y = "PC2") +
  theme(legend.position = "right")




######################PAIRWISE Comparison############################
library(pairwiseAdonis)


# Jaccard Distance Pairwise Comparisons
pairwise_jaccard <- pairwise.adonis(jaccard_dist, factors = village_factor, sim.function = vegdist, sim.method = "jaccard", p.adjust.m = "bonferroni")
print(pairwise_jaccard)

# Euclidean Distance Pairwise Comparisons
pairwise_euclidean <- pairwise.adonis( euclidean_dist, factors =  village_factor, sim.function = dist, sim.method = "euclidean", p.adjust.m = "bonferroni")
print(pairwise_euclidean)

# Manhattan Distance Pairwise Comparisons
pairwise_manhattan <- pairwise.adonis(manhattan_dist, factors = village_factor, sim.function = dist, sim.method = "manhattan", p.adjust.m = "bonferroni")
print(pairwise_manhattan)

