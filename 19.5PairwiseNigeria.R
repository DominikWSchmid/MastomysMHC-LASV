# Install and load pairwiseAdonis package if not already installed
if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
}
library(pairwiseAdonis)
library(reshape2)
library(ggplot2)

load("19.5PairwiseNigeria.RData")
# Jaccard Distance Pairwise Comparisons
pairwise_jaccard <- pairwise.adonis(jaccard_dist, factors = village_factor, sim.function = vegdist, sim.method = "jaccard", p.adjust.m = "bonferroni")
print(pairwise_jaccard)

# Euclidean Distance Pairwise Comparisons
pairwise_euclidean <- pairwise.adonis( euclidean_dist, factors =  village_factor, sim.function = dist, sim.method = "euclidean", p.adjust.m = "bonferroni")
print(pairwise_euclidean)

# Manhattan Distance Pairwise Comparisons
pairwise_manhattan <- pairwise.adonis(manhattan_dist, factors = village_factor, sim.function = dist, sim.method = "manhattan", p.adjust.m = "bonferroni")
print(pairwise_manhattan)


# Prepare matrix for heatmap


# Convert pairwise results to a matrix
pairwise_matrix <- matrix(NA, ncol = length(levels(village_factor)), nrow = length(levels(village_factor)),
                          dimnames = list(levels(village_factor), levels(village_factor)))



# Parse the village pairs and populate the matrix
for (i in 1:nrow(pairwise_jaccard)) {
  # Split the pair names (assuming they are stored in a single column)
  pair <- unlist(strsplit(as.character(pairwise_jaccard$pairs[i]), " vs "))
  
  # Assign the R² value to the corresponding matrix cells
  pairwise_matrix[pair[1], pair[2]] <- pairwise_jaccard$R2[i]
  pairwise_matrix[pair[2], pair[1]] <- pairwise_jaccard$R2[i] # Symmetric assignment
}

# View the matrix
pairwise_matrix

# Fill diagonals with 0 or NA if needed
diag(pairwise_matrix) <- NA

# Create heatmap
pairwise_df <- melt(pairwise_matrix, na.rm = TRUE)
ggplot(pairwise_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey") +
  labs(title = "Pairwise R² Heatmap", x = "Village", y = "Village", fill = "R²") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######################
# Create an empty list matrix to hold R² values and significance
pairwise_matrix <- matrix(list(), ncol = length(levels(village_factor)), nrow = length(levels(village_factor)),
                          dimnames = list(levels(village_factor), levels(village_factor)))

# Parse the village pairs and populate the matrix
for (i in 1:nrow(pairwise_jaccard)) {
  # Split the pair names (assuming they are stored in a single column)
  pair <- unlist(strsplit(as.character(pairwise_jaccard$pairs[i]), " vs "))
  
  # Create a string combining R² and significance
  result <- paste0("R²=", round(pairwise_jaccard$R2[i], 3), 
                   ", p=", pairwise_jaccard$p.adjusted[i], 
                   ifelse(pairwise_jaccard$sig[i] != "", paste0(" (", pairwise_jaccard$sig[i], ")"), ""))
  
  # Assign the result to the corresponding matrix cells
  pairwise_matrix[pair[1], pair[2]] <- result
  pairwise_matrix[pair[2], pair[1]] <- result # Symmetric assignment
}

# View the matrix
pairwise_matrix


#################


# Extract R² values as numeric for fill
r2_matrix <- matrix(as.numeric(gsub("^R²=([0-9\\.]+).*", "\\1", pairwise_matrix)), 
                    ncol = ncol(pairwise_matrix), 
                    nrow = nrow(pairwise_matrix),
                    dimnames = dimnames(pairwise_matrix))

# Create a data frame for ggplot
pairwise_df <- as.data.frame(as.table(r2_matrix))
colnames(pairwise_df) <- c("Village1", "Village2", "R2")


# Initialize matrix with NA_character_ instead of list()
pairwise_matrix <- matrix(NA_character_, 
                          ncol = length(levels(village_factor)), 
                          nrow = length(levels(village_factor)),
                          dimnames = list(levels(village_factor), levels(village_factor)))

# Fill diagonals with NA
diag(pairwise_matrix) <- NA_character_


for (i in 1:nrow(pairwise_jaccard)) {
  pair <- unlist(strsplit(as.character(pairwise_jaccard$pairs[i]), " vs "))
  
  result <- paste0("R²=", round(pairwise_jaccard$R2[i], 3), 
                   ", p=", pairwise_jaccard$p.adjusted[i], 
                   ifelse(pairwise_jaccard$sig[i] != "", paste0(" (", pairwise_jaccard$sig[i], ")"), ""))
  
  pairwise_matrix[pair[1], pair[2]] <- result
  pairwise_matrix[pair[2], pair[1]] <- result
}


pairwise_df$Label <- as.vector(pairwise_matrix)




# Plot heatmap


ggplot(pairwise_df, aes(x = Village1, y = Village2, fill = R2)) +
  geom_tile(color = "white") +                            # Heatmap tiles
  geom_text(aes(label = Label), color = "black", size = 3) + # Overlay labels
  scale_fill_gradient(low = "white", high = "red", na.value = "grey") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Pairwise R² Heatmap with Significance", 
       x = "Village 1", 
       y = "Village 2", 
       fill = "R² Value")


