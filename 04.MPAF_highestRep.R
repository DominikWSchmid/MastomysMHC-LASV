setwd("//134.60.87.178/Student Theses/Mastomys_Ayo/MHCI_2024_JanSarapak/R-Projects/2024_11_12_Post-Acacia-Final")

library(plyr)
library(dplyr)
library(readr)

library(ggplot2)


rm(list=ls())


c <-read.csv("Nigeria_Guinea_highestRep_XL.csv", header = TRUE)

  
c$ID <- as.factor(c$ID)  
  
  
  
c$SUM <- ave(c$COUNT, c$ID, FUN=sum)
c$RPAF <- c$COUNT/c$SUM #relative proportion of reads retrived for one particular allele 
#out of the total amount of reads returned within one given amplicon (relative per amplicon frequency)


c$ID <- as.factor(c$ID)



nID<-nlevels(c$ID)

c.MPAF.dat<-ddply(c, c("ALLELE"), summarise, meanMPAF=mean(RPAF), 
                  meanRead=mean(COUNT),
                  sdMPAF=sd(RPAF), 
                  n=length(ALLELE), 
                  se=sdMPAF/sqrt(n), 
                  prop=n/nID)  #mean per amplicon frequency, # prop n/number of ids
# number of IDs: 

# Singletons Drop rows where n == 1
c.MPAF.dat <- subset(c.MPAF.dat, n > 1)

# Calculate the total number of observations after filtering
total_obs <- sum(c.MPAF.dat$n)



c.reordered <- c.MPAF.dat[order(-c.MPAF.dat$meanMPAF),]



### 1. comparison of the two most common alleles with the remaining (following Huchard et al. 2012)###

c.reordered[,"rank"] <- c(1:407)



write_csv(c.reordered, "Nigeria_Guinea_Rep_reordered.csv")


str(c.reordered)

#We expected MPAF to be much higher for the two most common allels than for the remaining, 
#which would allow us to derive an estimate of the threshold MPAF 
#under which alleles are likely to represent artifactual alleles

# Arrange the data by 'prop' in descending order
c.reordered <- c.reordered %>%
  arrange(desc(prop))

# Select the second row and extract 'meanMPAF'
MPAF_Cutoff <- round(c.reordered$meanMPAF[2],3)

# Display the result
MPAF_Cutoff
# --> this is the MPAF cutoff 
# proportion cutoff is 0.01 




ggplot(c.reordered, aes(meanMPAF, prop, label=ALLELE))+
  geom_point(shape=21, size=1.5)+
  geom_point(data = c.reordered %>% filter(ALLELE == "ManaMHCI*001"), shape = 21, fill = "green") +
  geom_point(data = c.reordered %>% filter(ALLELE == "ManaMHCI*002"), shape = 21, fill = "green") +
  theme_bw()+
  geom_vline(xintercept = MPAF_Cutoff, linetype="dotted", color="red") +
  geom_hline(yintercept = 0.01, linetype="dotted", color="red") +
  xlim(0,0.4)+ylim(0,0.9)+
  ylab("Proportion")+xlab("Mean per amplicon frequency")+
  ggtitle("All_OnlyHighestReps")





######################## Apply Cutoff to whole Dataset


c <-read.csv("GuineaNigeria_allelereport_XL.csv", header = TRUE)

c$ID <- as.factor(c$ID)  



c$SUM <- ave(c$COUNT, c$ID, FUN=sum)
c$RPAF <- c$COUNT/c$SUM #relative proportion of reads retrived for one particular allele 
#out of the total amount of reads returned within one given amplicon (relative per amplicon frequency)

str(c)
c$ID <- as.factor(c$ID)

str(c)


nID<-nlevels(c$ID)

c.MPAF.dat<-ddply(c, c("ALLELE"), summarise, meanMPAF=mean(RPAF), 
                  meanRead=mean(COUNT),
                  sdMPAF=sd(RPAF), 
                  n=length(ALLELE), 
                  se=sdMPAF/sqrt(n), 
                  prop=n/nID)  #mean per amplicon frequency, # prop n/number of ids
# number of IDs: Guinea: 389, Nigeria: 589 , combined: 978

# Singletons Drop rows where n == 1
c.MPAF.dat <- subset(c.MPAF.dat, n > 1)

# Calculate the total number of observations after filtering
total_obs <- sum(c.MPAF.dat$n)



c.reordered <- c.MPAF.dat[order(-c.MPAF.dat$meanMPAF),]




### 1. comparison of the two most common alleles with the remaining (following Huchard et al. 2012)###

c.reordered[,"rank"] <- c(1:425)


write_csv(c.reordered, "Nigeria_Guinea_reordered.csv")





#We expected MPAF to be much higher for the two most common allels than for the remaining, 
#which would allow us to derive an estimate of the threshold MPAF 
#under which alleles are likely to represent artifactual alleles



ggplot(c.reordered, aes(meanMPAF, prop, label=ALLELE))+
  geom_point(shape=21, size=1.5)+
  geom_point(data = c.reordered %>% filter(ALLELE == "ManaMHCI*001"), shape = 21, fill = "green") +
  geom_point(data = c.reordered %>% filter(ALLELE == "ManaMHCI*002"), shape = 21, fill = "green") +
  theme_bw()+
  geom_vline(xintercept = MPAF_Cutoff, linetype="dotted", color="red") +
  geom_hline(yintercept = 0.01, linetype="dotted", color="red") +
  xlim(0,0.4)+ylim(0,0.9)+
  ylab("Proportion")+xlab("Mean per amplicon frequency")+
  ggtitle("All")

c.filtered<- c.reordered %>% 
  filter(meanMPAF >= MPAF_Cutoff | prop > 0.01)


ggplot(c.filtered, aes(meanMPAF, prop, label=ALLELE))+
  geom_point(shape=21, size=1.5)+
  geom_point(data = c.filtered %>% filter(ALLELE == "ManaMHCI*001"), shape = 21, fill = "green") +
  geom_point(data = c.filtered %>% filter(ALLELE == "ManaMHCI*002"), shape = 21, fill = "green") +
  theme_bw()+
  geom_vline(xintercept = MPAF_Cutoff, linetype="dotted", color="red") +
  geom_hline(yintercept = 0.01, linetype="dotted", color="red") +
  xlim(0,0.4)+ylim(0,0.9)+
  ylab("Proportion")+xlab("Mean per amplicon frequency")+
  ggtitle("All, filtered meanMPAF > 0.05 | prop > 0.01")



write_csv(c.filtered, "Nigeria_Guinea_filtered.csv")



# Step 1: Get the unique alleles in the original dataframe
original_alleles <- unique(c.reordered$ALLELE)

# Step 2: Get the unique alleles in the filtered dataframe
filtered_alleles <- unique(c.filtered$ALLELE)

# Step 3: Find the alleles that were removed (present in the original but not in the filtered dataframe)
alleles_removed <- setdiff(original_alleles, filtered_alleles)

# Print the removed alleles
cat("Alleles removed after filtering MPAF 0.045, p 0.01:", alleles_removed, "\n")



alleles_removed_Reps <- data.frame(allele=alleles_removed)

alleles_removed_all <- data.frame(allele = alleles_removed)

# Finding common alleles
common_alleles <- intersect(alleles_removed_Reps$allele, alleles_removed_all$allele)
common_alleles <- data.frame(allele = common_alleles)

# Finding unique alleles in each data frame
unique_to_Reps <- setdiff(alleles_removed_Reps$allele, alleles_removed_all$allele)
unique_to_all <- setdiff(alleles_removed_all$allele, alleles_removed_Reps$allele)

# Converting to data frames
unique_to_Reps <- data.frame(allele = unique_to_Reps)
unique_to_all <- data.frame(allele = unique_to_all)









b <-read.csv("GuineaNigeria_allelereport_XL.csv", header = TRUE)

b$ID <- as.factor(b$ID)  


# Step 1: Remove rows where the ALLELE is in the list of alleles_removed
b <- b %>% 
  filter(!ALLELE %in% alleles_removed)

# Step 2: Check the structure of the filtered dataframe
write_csv(b, "Nigeria_Guinea_XL_filtered.csv")

