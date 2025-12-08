
setwd("C:/Users/jansa/Desktop/2025-Mastomys-Final-Really")
rm(list=ls())

##### libraries for this script #####
library(cooccur)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lmerTest)
library(lme4)
library(hrbrthemes)
library(RColorBrewer)
library(tidyverse)

setwd("C:/Users/Jan_S/Desktop/2025-01-Mastomys-Final-Really")

load("21.Cooccurence.RData")

M <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv", header = TRUE)


threshold <- 0.1*nrow(M) 

# Identify ManaMHCI_ columns
mana_cols <- grepl("^ManaMHC_", colnames(M))

# Calculate column sums for these columns
col_sums <- colSums(M[, mana_cols])

# Select columns with a sum greater than the threshold
keep_cols <- names(col_sums[col_sums > threshold])

# Create a new dataset with the filtered columns
M <- M[, c(setdiff(colnames(M), names(col_sums)), keep_cols)]

str(M)
countries <- list(
  "Nigeria" = c("Ekpoma", "Ebudin", "Owo"),
  "Guinea" = c("Damania", "Sokourala", "Sonkonia")
)

marker_types <- list(
  "ST" = "^Supertype_",
  "AL" = "^ManaMHC_"
)

for (country in names(countries)) {
  # Subset metadata
  Mana <- subset(M, M$Village %in% countries[[country]])
  
  # Clean infection status
  Mana <- Mana %>%
    mutate(
      LASV_positive = ifelse(LASV_positive == "neg", 0, 1),
      IgG_positive = ifelse(IgG_positive == "neg", 0, 1),
      uninfected = ifelse(number_alleles == 1, 0, 1)
    )
  
  for (marker in names(marker_types)) {
    # Get marker columns
    marker_cols <- grep(marker_types[[marker]], colnames(Mana), value = TRUE)
    
    # Create co-occurrence matrix
    cooccur_Mana <- subset(Mana, select = c("LASV_positive", "IgG_positive", marker_cols))
    
    # Transpose matrix
    cooccur_Mana_long <- t(as.matrix(cooccur_Mana))
    cooccur_Mana_long <- cooccur_Mana_long[rowSums(cooccur_Mana_long) != 0, ]

    # Run co-occurrence
    ManaSTs <- cooccur(mat = cooccur_Mana_long, type = "spp_site", thresh = FALSE, spp_names = TRUE)
    
    
    # Create identifiers
    short_code <- paste0(substr(country, 1, 2), marker)  # e.g., "NiST"
    complete_name <- paste0("complete_", short_code)
    
    # Run full summary table
    complete_ManaSTs <- prob.table(ManaSTs)
    effectsize_Mana <- effect.sizes(ManaSTs, standardized = TRUE, matrix = FALSE)
    complete_ManaSTs$effects <- effectsize_Mana$effects
    complete_ManaSTs$obs.v.exp <- (complete_ManaSTs$exp_cooccur - complete_ManaSTs$obs_cooccur) * -1
    
    # Assign to global environment
    assign(short_code, ManaSTs, envir = .GlobalEnv)
    assign(complete_name, complete_ManaSTs, envir = .GlobalEnv)
    

  }
}
################################## Nigeria Supertypes  ################
summary(NiST)
plot(NiST)

pair(NiST, "LASV_positive")
pair(NiST, "IgG_positive")


print(NiST)
print(complete_NiST)
obs.v.exp(NiST) 


######### Nigeria Alleles
print(paste("N markers found for", short_code, ":", length(marker_cols)))


summary(NiAL)
plot(NiAL)

pair(NiAL, "LASV_positive")
pair(NiAL, "IgG_positive")


print(NiAL)
print(complete_NiAL)
obs.v.exp(NiAL) 


################################# Guinea Supertypes ############
summary(GuST)
plot(GuST)

pair(GuST, "LASV_positive")
pair(GuST, "IgG_positive")


print(GuST)
print(complete_GuST)
obs.v.exp(GuST) 


######### Guineea Alleles 
summary(GuAL)
plot(GuAL)

pair(GuAL, "LASV_positive")
pair(GuAL, "IgG_positive")


print(GuAL)
print(complete_GuAL)
obs.v.exp(GuAL) 

##########################  Allele Correlation for r > 0.08 #########

# Datenaufbereitung für jedes Land
mana_data_list <- list()
for (country in names(countries)) {
  # Subset metadata
  Mana <- subset(M, M$Village %in% countries[[country]])
  
  # Clean infection status
  Mana <- Mana %>%
    mutate(
      LASV_positive = ifelse(LASV_positive == "neg", 0, 1),
      IgG_positive = ifelse(IgG_positive == "neg", 0, 1),
      uninfected = ifelse(number_alleles == 1, 0, 1)
    )
  mana_data_list[[country]] <- Mana
}

# Separate Korrelationsanalyse der Allele
allele_correlation_results <- list()
for (country in names(mana_data_list)) {
  Mana <- mana_data_list[[country]]
  allele_data <- Mana %>% select(starts_with("ManaMHC_"))
  
  if (ncol(allele_data) > 1) {
    allele_cor_matrix <- cor(allele_data)
    high_cor_alleles <- which(abs(allele_cor_matrix) >= 0.8 & upper.tri(allele_cor_matrix), arr.ind = TRUE)
    
    if (nrow(high_cor_alleles) > 0) {
      high_cor_pairs <- data.frame(
        Allele1 = colnames(allele_data)[high_cor_alleles[, 2]],
        Allele2 = colnames(allele_data)[high_cor_alleles[, 1]],
        Correlation = allele_cor_matrix[high_cor_alleles]
      )
      
      allele_frequencies <- colSums(allele_data)
      most_frequent_alleles <- apply(high_cor_pairs[, c("Allele1", "Allele2")], 1, function(pair) {
        freq1 <- allele_frequencies[pair[1]]
        freq2 <- allele_frequencies[pair[2]]
        if (freq1 >= freq2) return(pair[1]) else return(pair[2])
      })
      high_cor_pairs$MostFrequent <- most_frequent_alleles
      allele_correlation_results[[country]] <- high_cor_pairs
      print(paste("Hochkorrelierte Allelpaare in", country, ":"))
      print(high_cor_pairs)
    } else {
      print(paste("Keine hochkorrelierten Allelpaare (r >= 0.8) in", country))
      allele_correlation_results[[country]] <- data.frame() # Leerer Dataframe, falls keine Korrelationen
    }
  } else {
    print(paste("Nicht genügend Allele für die Korrelationsanalyse in", country))
    allele_correlation_results[[country]] <- data.frame()
  }
}

# Die ursprüngliche Co-Okkurrenzanalyse bleibt in der zweiten Schleife
for (country in names(countries)) {
  # Subset metadata (verwenden jetzt die vorbereiteten Daten)
  Mana <- mana_data_list[[country]]
  
  for (marker in names(marker_types)) {
    # Get marker columns
    marker_cols <- grep(marker_types[[marker]], colnames(Mana), value = TRUE)
    
    # Create co-occurrence matrix
    cooccur_Mana <- subset(Mana, select = c("LASV_positive", "IgG_positive", marker_cols))
    
    # Transpose matrix for co-occurrence
    cooccur_Mana_long <- t(as.matrix(cooccur_Mana))
    cooccur_Mana_long <- cooccur_Mana_long[rowSums(cooccur_Mana_long) != 0, ]
    
    # Run co-occurrence
    ManaSTs <- cooccur(mat = cooccur_Mana_long, type = "spp_site", thresh = FALSE, spp_names = TRUE)
    
    # Create identifiers
    short_code <- paste0(substr(country, 1, 2), marker)  # e.g., "NiST"
    complete_name <- paste0("complete_", short_code)
    
    # Run full summary table
    complete_ManaSTs <- prob.table(ManaSTs)
    effectsize_Mana <- effect.sizes(ManaSTs, standardized = TRUE, matrix = FALSE)
    complete_ManaSTs$effects <- effectsize_Mana$effects
    complete_ManaSTs$obs.v.exp <- (complete_ManaSTs$exp_cooccur - complete_ManaSTs$obs_cooccur) * -1
    
    # Assign co-occurrence results to global environment
    assign(short_code, ManaSTs, envir = .GlobalEnv)
    assign(complete_name, complete_ManaSTs, envir = .GlobalEnv)
  }
}

print(allele_correlation_results)


######################### HEATMAPS #############


heatmap.NiST<-subset(complete_NiST, sp1_name=="LASV_positive" | sp1_name=="IgG_positive"  
                     #| sp1_name=="uninfected"
)
heatmap.NiST<-subset(heatmap.NiST, sp2_name!="IgG_positive")


heatmap.GuST<-subset(complete_GuST, sp1_name=="LASV_positive" | sp1_name=="IgG_positive"  
                    #| sp1_name=="uninfected"
)
heatmap.GuST<-subset(heatmap.GuST, sp2_name!="IgG_positive")


heatmap.NiST <- heatmap.NiST %>%
  mutate(Asterisk = ifelse(p_lt <= 0.01 | p_gt <= 0.01, "*", NA),
         Country = "Nigeria")

heatmap.GuST <- heatmap.GuST %>%
  mutate(Asterisk = ifelse(p_lt <= 0.01 | p_gt <= 0.01, "*", NA),
         Country = "Guinea")

################ LOOP  ###########


# Define your datasets and corresponding labels
heatmap_sets <- list(
  NiST = complete_NiST,
  GuST = complete_GuST,
  NiAL = complete_NiAL,
  GuAL = complete_GuAL
)

# Prepare an empty list to hold processed dfs
heatmap_results <- list()

for (name in names(heatmap_sets)) {
  df <- heatmap_sets[[name]]
  
  # Filter for relevant infection status
  df_filtered <- df %>%
    filter(sp1_name %in% c("LASV_positive", "IgG_positive"),
           sp2_name != "IgG_positive") %>%
    mutate(
      Asterisk = ifelse(p_lt <= 0.01 | p_gt <= 0.01, "*", NA),
      Country = ifelse(grepl("^Ni", name), "Nigeria", "Guinea"),
      Marker = substr(name, 3, 4)  # "ST" or "AL"
    )
  
  heatmap_results[[name]] <- df_filtered
}






# Combine ST (Supertypes) datasets
heatmap_combined_ST <- rbind(
  heatmap_results$NiST,
  heatmap_results$GuST
)

# Combine AL (Alleles) datasets
heatmap_combined_AL <- rbind(
  heatmap_results$NiAL,
  heatmap_results$GuAL
)


########## Heatmap sT


# Find combinations where both countries show significance
shared_significance_ST <- heatmap_combined_ST %>%
  filter(!is.na(Asterisk)) %>%
  group_by(sp1_name, sp2_name) %>%
  tally() %>%
  filter(n == 2) %>%
  mutate(shared_sig = TRUE)

# Merge back to flag tiles
heatmap_combined_ST <- heatmap_combined_ST %>%
  left_join(shared_significance_ST, by = c("sp1_name", "sp2_name")) %>%
  mutate(shared_sig = ifelse(is.na(shared_sig), FALSE, shared_sig))

heatmap_combined_ST <- heatmap_combined_ST %>%
  mutate(
    border_color = ifelse(shared_sig, "black", "grey50"),
    border_width = ifelse(shared_sig, 1, 1)
  )

ggplot(heatmap_combined_ST, aes(sp2_name, sp1_name, fill = effects)) + 
  geom_tile(aes(color = border_color, linewidth = border_width)) +
  scale_color_identity() +  # uses the actual hex values/names from the column
  scale_linewidth_identity() +  # same for linewidth
  scale_fill_gradient2(high = "red", mid = "white", low = "darkgreen", 
                       limits = c(-0.05, 0.05)) +
    geom_text(aes(label = Asterisk), size = 7, vjust = .7, hjust = 0.5) +
  facet_wrap(~Country, ncol = 1) +
  scale_x_discrete(
    limits = paste0("Supertype_", 1:19),
    labels = paste0("ST", 1:19)  ) +  
  scale_y_discrete(    labels = c("LASV_positive" = "LASV+", "IgG_positive" = "IgG+")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14, face = "bold")) +
  ylab("Infection status") + xlab("Supertype") +
  coord_fixed(ratio = 1)+
  ggtitle("p<0.01 = *")


############ HEATMAP ALLELE


# Find combinations where both countries show significance
shared_significance_AL <- heatmap_combined_AL %>%
  filter(!is.na(Asterisk)) %>%
  group_by(sp1_name, sp2_name) %>%
  tally() %>%
  filter(n == 2) %>%
  mutate(shared_sig = TRUE)

# Merge back to flag tiles
heatmap_combined_AL <- heatmap_combined_AL %>%
  left_join(shared_significance_AL, by = c("sp1_name", "sp2_name")) %>%
  mutate(shared_sig = ifelse(is.na(shared_sig), FALSE, shared_sig))

heatmap_combined_AL <- heatmap_combined_AL %>%
  mutate(
    border_color = ifelse(shared_sig, "black", "grey50"),
    border_width = ifelse(shared_sig, 1, .5)
  )

heatmap_combined_AL <- heatmap_combined_AL %>%
  mutate(
    allele_number = gsub(".*_", "", sp2_name),             # keep only number part
    allele_number = as.numeric(allele_number)              # turn to numeric for proper sorting
  ) %>%
  arrange(allele_number) %>%
  mutate(
    sp2_name = sprintf("%03d", allele_number)              # overwrite to 3-digit style (e.g. 001, 019)
  )


ggplot(heatmap_combined_AL, aes(sp2_name, sp1_name, fill = effects)) + 
  geom_tile(aes(color = border_color, linewidth = border_width)) +
  scale_color_identity() + 
  scale_linewidth_identity() + 
  scale_fill_gradient2(
    high = "red", mid = "white", low = "darkgreen",     limits = c(-0.08, 0.08)
  ) +
  geom_text(aes(label = Asterisk), size = 7, vjust = .7, hjust = 0.5) +
  facet_wrap(~Country, ncol = 1) +
  
  # Optional: reorder or relabel alleles if needed
  # scale_x_discrete(limits = sort(unique(heatmap_combined_AL$sp2_name))) +
  
  scale_y_discrete(
    labels = c("LASV_positive" = "LASV+", "IgG_positive" = "IgG+")
  ) +
  scale_x_discrete()+
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  
  ylab("Infection status") + xlab("Allele") +
  coord_fixed(ratio = 1)#+
  #ggtitle("p<0.01 = *")


############## Only significant Allele

# Get alleles with significance in at least one country
significant_alleles <- heatmap_combined_AL %>%
  filter(!is.na(Asterisk)) %>%
  pull(sp2_name) %>%
  unique()

# Filter main dataframe
filtered_heatmap_AL <- heatmap_combined_AL %>%
  filter(sp2_name %in% significant_alleles)


heatmap<-ggplot(filtered_heatmap_AL, aes(sp2_name, sp1_name, fill = effects)) + 
  geom_tile(aes(color = border_color, linewidth = border_width)) +
  scale_color_identity() + 
  scale_linewidth_identity() + 
  scale_fill_gradient2(
    high = "red", mid = "white", low = "darkgreen",
    limits = c(-0.08, 0.08)
  ) +
  geom_text(aes(label = Asterisk), size = 7, vjust = .7, hjust = 0.5) +
  facet_wrap(~Country, ncol = 1) +
  scale_y_discrete(
    labels = c("LASV_positive" = "LASV+", "IgG_positive" = "IgG+")
  ) +
  scale_x_discrete(limits = sort(unique(filtered_heatmap_AL$sp2_name))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  ylab("Infection status") + xlab("Mana-DRB*") +
  coord_fixed(ratio = 1)+theme(axis.title = element_text(face="bold", size=14))
   
####### table output
allele_cooccurrence_table <- filtered_heatmap_AL %>%
  select(Allele = sp2_name,
         Infection_Status = sp1_name,
         Country,
         Effect_Size = effects,
         P_Value_Significant = Asterisk) %>%
  arrange(Country, Infection_Status, Allele)

print(allele_cooccurrence_table)

