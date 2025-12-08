setwd("//134.60.87.178/f/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")



library(dplyr)
library(readxl)
library(lubridate) 


MHC_OneSheet <- read_excel("MHC_OneSheet.xlsx", na = "NA")
MHC_OneSheet$Village[MHC_OneSheet$Village == "Eguare-Egoro"] <- "Ekpoma"

os <- MHC_OneSheet



wide_data<- read.csv("Nigeria_Guinea_wide_allele.csv")



wide_data <- wide_data %>%
  mutate(Run = ifelse(grepl("-j$", ID), "JS", Run))

wide_data$ID <- gsub("-j", "", wide_data$ID)
wide_data$ID <- gsub(("-rep|-rep2"), "", wide_data$ID)

# Assuming os$IgG_postive is your factor variable
# Step 1: Convert the factor to a character vector
os$igG_char <- as.character(os$IgG_positive)

# Step 2: Replace "." with NA
os$igG_char[os$igG_char == "."] <- NA



names(os)[names(os) == "Sample_No."] <- "ID"
# Example data frames
# os <- data.frame(ID = c(...), LASV_positive = c(...), IgG_positive = c(...))
# wide_data <- data.frame(ID = c(...), other_columns = c(...))
#### Add Sex, Wt, ELW to metadata

NG_S <- read_excel("AutopsiesDataNigeria.xlsx", na = "NA")
NG_S <- read_excel("SexWeightEyelenseweight_Nigeria.xlsx", na = "NA")

# Ensure column names align for merging
colnames(NG_S)[colnames(NG_S) == "N viro"] <- "ID"

# Select only relevant columns for merging
NG_S_subset <- NG_S %>%
  select(ID, ELW, Wt, Sex, HB)

# Merge the new data with the existing data frame
os <- merge(os, NG_S_subset, by = "ID", all.x = TRUE)




##### DATE 

os$Date <- gsub("Sept,", "Sep,", os$Date) # Fix "Sept" to "Sep"
os$Date <- gsub("\\s+", " ", os$Date)    # Remove extra spaces
os$Date <- gsub("\\.0$", "", os$Date) 

os$Date_Standardized <- sapply(os$Date, function(x) {
  if (grepl("^\\d+$", x)) {
    as.character(as.Date(as.numeric(x), origin = "1899-12-30"))
  } else {
    parsed <- parse_date_time(x, orders = c("dmy", "mdy", "ymd", "d B Y", "A, d B, Y"))
    if (is.na(parsed)) {
      NA  # Fallback if parsing still fails
    } else {
      as.character(parsed)
    }
  }
})

os$Country <- "Nigeria"


# Merge the data frames based on the ID column
merged_data <- merge(wide_data, os[, c("ID", "LASV_positive", "IgG_positive", "Country", "State", "Village", "Date_Standardized", "Sex", "Wt", "ELW", "HB")], by = "ID", all.x = TRUE)




# View the result to ensure columns are correctly reordered
head(merged_data)



coordinates <- read_excel("CoordinatesNigeriaGuinea.xlsx", na = "NA")
merged_df <- merge(merged_data,coordinates [, c("Village","lattitude","longitude")],by = "Village", all.x = T)

# Reorder columns: move LASV_positive and IgG_positive to second and third positions
merged_data <- merged_df[c("ID", "LASV_positive", "IgG_positive", "Country", "State", "Village","Date_Standardized" , "Sex", "Wt", "ELW", "HB", "lattitude", "longitude", 
                             setdiff(names(merged_data), c("ID", "LASV_positive", "IgG_positive", "Country", "State", "Village","Date_Standardized" , "Sex", "Wt", "ELW","HB", "lattitude", "longitude")))]

# write.csv(merged_data, "Nigeria_Guinea_Wide.csv", row.names = FALSE)

### GUINEAN




guinea_meta <- read_excel("MHC_Guinea.xlsx", na = "NA")

guinea_meta$ELW <- as.numeric(guinea_meta$ELW)
    
str(merged_data)

str(guinea_meta)

colnames(guinea_meta)[colnames(guinea_meta) == "Nviro"] <- "ID"

# Step 1: Ensure first three characters (prefix) are capital letters
guinea_meta$ID <- sub("^(\\w)(\\w+)", "\\U\\1\\L\\2", guinea_meta$ID, perl = TRUE)

# Step 2: Replace the blank space with an underscore
guinea_meta$ID <- gsub(" ", "-", guinea_meta$ID)

# Step 3: Ensure that the numbers are three digits
guinea_meta$ID <- gsub("(\\D)(\\d{1})(\\D|$)", "\\100\\2\\3", guinea_meta$ID)
guinea_meta$ID <- gsub("(\\D)(\\d{2})(\\D|$)", "\\10\\2\\3", guinea_meta$ID)

# Display the result
guinea_meta$ID

merged_data$ID



# Step 1: Identify and modify only IDs that have a prefix (e.g., Dam, SOK, TAM, SON)
merged_data$ID <- sapply(merged_data$ID, function(id) {
  # Check if ID contains a prefix like "Dam", "SOK", "TAM", etc.
  if (grepl("^[A-Za-z]+", id)) {
    
    # Capitalize first letter, lowercase the next two letters (prefix)
    id <- sub("^(\\w)(\\w+)", "\\U\\1\\L\\2", id, perl = TRUE)
    
    # Replace space or hyphen with an underscore
    id <- gsub("[ -]", "-", id)
    
    # Add leading zeros to make the number three digits, while keeping '-rep' part unchanged
    id <- gsub("(\\D)(\\d{1})(\\D|$)", "\\100\\2\\3", id)
    id <- gsub("(\\D)(\\d{2})(\\D|$)", "\\10\\2\\3", id)
  }
  
  # Step 2: Return the ID unchanged if it's numeric-only or has '-rep' suffix
  return(id)
})

# Display the modified IDs
merged_data$ID


str(guinea_meta)



####### Date 

guinea_meta$Date <- guinea_meta$date
guinea_meta$Date <- gsub("Sept,", "Sep,", guinea_meta$Date) # Fix "Sept" to "Sep"
guinea_meta$Date <- gsub("\\s+", " ", guinea_meta$Date)    # Remove extra spaces
guinea_meta$Date <- gsub("\\.0$", "", guinea_meta$Date) 

guinea_meta$Date_StandardizedG <- sapply(guinea_meta$Date, function(x) {
  if (grepl("^\\d+$", x)) {
    as.character(as.Date(as.numeric(x), origin = "1899-12-30"))
  } else {
    parsed <- parse_date_time(x, orders = c("dmy", "mdy", "ymd", "d B Y", "A, d B, Y"))
    if (is.na(parsed)) {
      NA  # Fallback if parsing still fails
    } else {
      as.character(parsed)
    }
  }
})

guinea_meta$CountryG <- "Guinea"


# Merge the two data frames by 'ID', keeping all rows from merged_d
merged_new <- merge(merged_data, guinea_meta[, c("ID", "LASV(ifa)", "LASV(pcr)", "CountryG", "Date_StandardizedG", "Wt", "Sex", "ELW", "HB")], by = "ID", all.x = TRUE)



# Update 'IgG_positive' only where it's NA, taking values from 'LASV(ifa)'
merged_new$IgG_positive <- ifelse(is.na(merged_new$IgG_positive), merged_new$`LASV(ifa)`, merged_new$IgG_positive)

# Update 'LASV_positive' only where it's NA, taking values from 'LASV(pcr)'
merged_new$LASV_positive <- ifelse(is.na(merged_new$LASV_positive), merged_new$`LASV(pcr)`, merged_new$LASV_positive)


# Update 'LASV_positive' only where it's NA, taking values from 'LASV(pcr)'
merged_new$Country <- ifelse(is.na(merged_new$Country), merged_new$`CountryG`, merged_new$Country)

merged_new$Date_Standardized <- ifelse(is.na(merged_new$Date_Standardized), merged_new$`Date_StandardizedG`, merged_new$Date_Standardized)


colnames(merged_new)[colnames(merged_new) == "Date_Standardized"] <- "Date_Capture"


# Remove the now unnecessary columns 'LASV(ifa)' and 'LASV(pcr)'
merged_new <- merged_new[, !colnames(merged_new) %in% c("LASV(ifa)", "LASV(pcr)", "CountryG", "Date_StandardizedG")]

# Display the first few rows of the new data frame
head(merged_new)
# view(merged_new)


# Use mutate and case_when to update the Village column based on the ID prefixes
merged_new <- merged_new %>%
  mutate(Village = case_when(
    grepl("^Dam", ID) ~ "Damania",
    grepl("^Tam", ID) ~ "Tambaya",
    grepl("^Son", ID) ~ "Sonkonia",
    grepl("^Sok", ID) ~ "Sokourala",
    TRUE ~ Village  # Keep the original Village value if no match
  ))
# Display the updated data frame
#View(merged_new)


merged_new <- merged_new %>%
  mutate(State = if_else(Run == "GN", "Faranah", State))



# Merge the two data frames
merged_with_coords <- merged_new %>%
  left_join(coordinates, by = "Village", suffix = c("", ".new")) %>%
  mutate(
    # Use new latitude and longitude if the original values are NA
    lattitude = coalesce(lattitude, lattitude.new),
    longitude = coalesce(longitude, longitude.new)
  ) %>%
  # Remove the new latitude and longitude columns after merging
  select(-lattitude.new, -longitude.new)



merged_with_coords$Run <- gsub("GN", "Guinea", merged_with_coords$Run)
merged_with_coords$Run <- gsub("NG", "Nigeria", merged_with_coords$Run)
merged_with_coords$Run <- gsub("JS", "Nigeria", merged_with_coords$Run)

############################### Sex Wt ELW

merged_with_coords <- merged_with_coords %>%
  mutate(
    Wt.x = as.numeric(Wt.x), # Wenn Wt.x der Character-Typ ist
    Wt.y = as.numeric(Wt.y), # Wenn Wt.y der Character-Typ ist
    Sex.x = as.character(Sex.x), # Stellen Sie sicher, dass beide 'Sex'-Spalten Character sind
    Sex.y = as.character(Sex.y),
    ELW.x = as.numeric(ELW.x), # Wenn ELW.x der Character-Typ ist
    ELW.y = as.numeric(ELW.y),  # Wenn ELW.y der Character-Typ ist
    HB.x = as.numeric(HB.x),
    HB.y = as.numeric(HB.y)
  ) %>%
  mutate(
    Wt = coalesce(Wt.x, Wt.y),
    Sex = coalesce(Sex.x, Sex.y),
    ELW = coalesce(ELW.x, ELW.y),
    HB = coalesce(HB.x, HB.y)
  ) %>%
  select(-Wt.x, -Wt.y, -Sex.x, -Sex.y, -ELW.x, -ELW.y, -HB.x, -HB.y)

merged_with_coords <- merged_with_coords %>%
  select(ID, LASV_positive, IgG_positive, Country, State, Village, Date_Capture, 
         Wt, Sex, ELW, HB, everything())  # Moves Wt, Sex, and ELW after Date_Capture


write.csv(merged_with_coords, "Nigeria_Guinea_Wide.csv", row.names = FALSE)

# Convert columns to factors
merged_with_coords$Acacia_Run <- as.factor(merged_with_coords$Run)
merged_with_coords$Country <- as.factor(merged_with_coords$Country)
merged_with_coords$Village <- as.factor(merged_with_coords$Village)


# Count values for 'Run'
run_counts <- table(merged_with_coords$AcaciaRun)

# Count values for 'Country'
country_counts <- table(merged_with_coords$Country)

# Count values for 'Village'
village_counts <- table(merged_with_coords$Village)

# Display the counts
run_counts
country_counts
village_counts



