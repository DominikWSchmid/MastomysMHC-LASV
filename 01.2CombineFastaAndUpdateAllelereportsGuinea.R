setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

library (BiocManager)
library (Biostrings)


# Replace 'first_run.fasta' and 'second_run.fasta' with your actual file names
first_run <- readDNAStringSet("Nigeriaupdate.fasta")
second_run <- readDNAStringSet("Guinea.fasta")

# Initialize a counter for the new sequence names
last_number <- length(first_run)
 
# Create a vector to keep track of name changes
name_changes <- c()

# Function to generate new names based on a counter
generate_new_name <- function(counter) {
  paste0("ManaMHCI*", counter)
}

# Create a new DNAStringSet for the combined sequences
combined_sequences <- first_run

# Match and rename sequences
for (i in seq_along(second_run)) {
  seq <- second_run[i]
  match_idx <- which(as.character(first_run) == as.character(seq))
  
  if (length(match_idx) > 0) {
    # If a match is found, rename the sequence to the corresponding name in the first run
    new_name <- names(first_run)[match_idx]
    names(seq) <- new_name
    name_changes <- c(name_changes, paste(names(second_run)[i], "->", new_name))
  } else {
    # If no match is found, generate a new name and add to the combined set
    last_number <- last_number + 1
    new_name <- generate_new_name(last_number)
    names(seq) <- new_name
    name_changes <- c(name_changes, paste(names(second_run)[i], "->", new_name))
    combined_sequences <- c(combined_sequences, seq)
  }
}

# Print the name changes
print(name_changes)

# Save the combined sequences to a new FASTA file
writeXStringSet(combined_sequences, "GuineaNigeria.fasta")




########################### CHANGE ALLELE NAMES IN Allelereports



library(readr)
library(dplyr)

# Create a named vector for easy lookup from name_changes
name_changes_vector <- setNames(
  sapply(strsplit(name_changes, " -> "), `[`, 2),
  sapply(strsplit(name_changes, " -> "), `[`, 1)
)



# Function to find the start of the second part of the CSV
find_second_part_start <- function(lines) {
  for (i in seq_along(lines)) {
    if (grepl("^Individual Report:", lines[i])) {
      return(i)
    }
  }
  return(NA)
}

# Read the entire CSV file of the new run as text

Ayonew_allelereport <- read_csv("Guinea_allelereport.csv", 
                                skip = 3)

file_lines <- readLines("Guinea_allelereport.csv")

# Find where the second part starts
second_part_start <- find_second_part_start(file_lines)


# Split the file into two parts if the second part exists
if (!is.na(second_part_start)) {
  first_part_lines <- file_lines[1:(second_part_start - 1)]
  second_part_lines <- file_lines[second_part_start:length(file_lines)]
  
  # Read the first part into a dataframe
  first_part <- read_csv(I(paste(first_part_lines, collapse = "\n")), col_names = c("Allele", "Abundance", "Sequence"), skip = 4)
  
  # Read the second part into a dataframe
  second_part <- read_csv(I(paste(second_part_lines, collapse = "\n")), col_names = c("Individual", "Nr_of_alleles", "Alleles"), skip = 2)
  
  # Update the Allele column in the first part based on the name_changes vector
  first_part$Allele <- ifelse(first_part$Allele %in% names(name_changes_vector),
                              name_changes_vector[first_part$Allele],
                              first_part$Allele)
  
  # Function to rename alleles in the second part
  rename_alleles <- function(allele_string, name_changes_vector) {
    alleles <- unlist(strsplit(allele_string, ";"))
    renamed_alleles <- sapply(alleles, function(x) ifelse(x %in% names(name_changes_vector), name_changes_vector[x], x))
    return(paste(renamed_alleles, collapse = ";"))
  }
  
  # Update the Alleles column in the second part based on the name_changes vector
  second_part$Alleles <- sapply(second_part$Alleles, rename_alleles, name_changes_vector = name_changes_vector)
  
 
  
  
  # Combine the first part and the second part back into one dataframe-like format
  combined_data <- c(
    first_part_lines[1:4], # Including headers of the first part
    paste(first_part$Allele, first_part$Abundance, first_part$Sequence, sep = ","),
    "", # Add an empty line to separate parts
    second_part_lines[1:2], # Including headers of the second part
    paste(second_part$Individual, second_part$Nr_of_alleles, second_part$Alleles, sep = ",")
  )
  
} else {
  stop("The second part 'Individual Report' was not found in the file.")
}



combined_dataf <- as.data.frame(combined_data)

combined_dataf[203,]

# write the combined data to a CSV file
writeLines(combined_data, "updated_Guinea_allelereport.csv")

####################### ALLELEREPORT_XL###################

# Read the new CSV file where sequence names need to be updated
allele_report_XL <- read_csv("Guinea_allelereport_XL.csv")

# Update the ALLELE column based on the name_changes vector
allele_report_XL <- allele_report_XL %>%
  mutate(ALLELE = ifelse(ALLELE %in% names(name_changes_vector), name_changes_vector[ALLELE], ALLELE))

allele_report_XL$ID <- gsub("SOK-02-rep2", "SOK-02-rep-j", allele_report_XL$ID, perl = TRUE)


# Function to normalize ID numbering to three digits
normalize_id <- function(id) {
  match <- regexpr("(\\d+)", id)  # Find the numeric part
  if (match > 0) {
    number <- as.numeric(regmatches(id, match))  # Extract and convert to numeric
    formatted_number <- sprintf("%03d", number)  # Format as three-digit
    id <- sub("(\\d+)", formatted_number, id)  # Replace in original string
  }
  return(id)
}

# Apply function to all IDs
allele_report_XL$ID <- sapply(allele_report_XL$ID, normalize_id)




# Write the updated dataframe back to a CSV file
write_csv(allele_report_XL, "updated_Guinea_allelereport_XL.csv")

