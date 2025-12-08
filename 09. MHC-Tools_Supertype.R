setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")

library(seqinr)
library(MHCtools)

library (BiocManager)
library (Biostrings)

library(dada2)






####### DISTCALC ()


dist_out <- "//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final/SupertypeOutput"
seq_file <-  read.fasta("GuineaNigeria_filtered_noGaps.fasta")




DistCalc(
  seq_file,
  path_out = dist_out,
  input_fasta = T,
  input_seq = "nucl",
  aa_dist = T,
  codon_pos = c(3, 6, 18, 39, 44, 56, 57, 59, 60, 63, 64, 67, 71, 77),
  dist_type = "S"
)



################# BOOTKMEANS()



# Path to your output matrix file
dist_matrix_path <- "2024_11_12_Post-Acacia-Final/SupertypeOutput/distance_matrix.csv"

list.files(path = path_out, full.names = TRUE)



# Paths to each file
z1_matrix_file <- "SupertypeOutput/z1_matrix_20241113.csv"
z2_matrix_file <- "SupertypeOutput/z2_matrix_20241113.csv"
z3_matrix_file <- "SupertypeOutput/z3_matrix_20241113.csv"
z4_matrix_file <- "SupertypeOutput/z4_matrix_20241113.csv"
z5_matrix_file <- "SupertypeOutput/z5_matrix_20241113.csv"

# Load the z matrices
z1_matrix <- as.matrix(read.csv(z1_matrix_file, row.names = 1))
z2_matrix <- as.matrix(read.csv(z2_matrix_file, row.names = 1))
z3_matrix <- as.matrix(read.csv(z3_matrix_file, row.names = 1))
z4_matrix <- as.matrix(read.csv(z4_matrix_file, row.names = 1))
z5_matrix <- as.matrix(read.csv(z5_matrix_file, row.names = 1))

# Output path for BootKmeans results
path_out <- "//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final/SupertypeOutput/BootKmeansResult"




# Run BootKmeans 
BootKmeans(
  z1_matrix, z2_matrix, z3_matrix, z4_matrix, z5_matrix,
  threshold = 0.01,
  no_scans = 1000,
  max_k = 40,
  iter.max = 1000000,
  nstart = 200,
  algorithm = "Hartigan-Wong",
  path_out = path_out
)








######### CLUSTERMATCH()


k_summary_table <- read.csv("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final/SupertypeOutput/BootKmeansResult/k_means_bootstrap_summary_stats_20241113.csv",
                            row.names = 1)

filepath  <- "//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final/SupertypeOutput/BootKmeansResult/Clusters"

clusterMatch_out <- "//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final/SupertypeOutput/Clustermatch"


ClusterMatch(filepath = filepath,
             path_out = clusterMatch_out, 
             k_summary_table = k_summary_table)



clusters <- read.csv("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final/SupertypeOutput/Clustermatch/ClusterMatch_summary_stats_20241114.csv")


clusters
