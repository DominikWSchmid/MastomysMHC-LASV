setwd("C:/Users/Jan_S/Desktop/2024_11_12_Post-Acacia-Final")



# Load required libraries
library(ape)
library(ggtree)
library(dplyr)
library(Biostrings)
library(ggplot2)
library(ggtreeExtra)

# Step 1: Read the Newick tree
tree <- read.tree("GuineaNigeria_filtered_tree.newick")

# Step 2: Read the FASTA file and extract allele names
fasta <- readDNAStringSet("GuineaNigeria_filtered_noGaps.fasta")
alleles_in_fasta <- names(fasta)


# Read the metadata
df <- read.csv("Nigeria_Guinea_Wide_ST.csv")


# Step 3: Calculate occurrence of alleles in each country
allele_occurrence <- df %>%
  select(Country, starts_with("ManaMHCI")) %>%
  pivot_longer(-Country, names_to = "Allele", values_to = "Present") %>%
  filter(Present > 0) %>%
  count(Allele, Country) %>%
  pivot_wider(names_from = Country, values_from = n, values_fill = 0)

# Step 4: Match occurrence data to tree tips
tree_data <- data.frame(Allele = tree$tip.label) %>%
  left_join(allele_occurrence, by = "Allele")

# Step 5: Visualize the tree with occurrence as a heatmap
# Visualize the tree with bar plots for occurrence
p <- ggtree(tree) +
  geom_tiplab(size = 2, align = TRUE) +  # Tip labels
  geom_fruit(
    data = tree_data,
    geom = geom_bar,
    aes(y = Allele, x = Guinea, fill = "Guinea"),  # Map 'y' explicitly
    width = 0.7,
    stat = "identity"
  ) +
  geom_fruit(
    data = tree_data,
    geom = geom_bar,
    aes(y = Allele, x = Nigeria, fill = "Nigeria"),  # Map 'y' explicitly
    width = 0.7,
    stat = "identity"
  ) +
  scale_fill_manual(values = c("Nigeria" = "#48ddd1", "Guinea" = "#ff725d")) +
  theme(legend.position = "right") +
  labs(title = "Phylogenetic Tree with Allele Occurrence by Country")

print(p)


