################# Ghana porject #########################
### cooccurrence analysis ###
### investigator: Dominik & Magdalena
### initial R script build by Tatiana & Thomas

#####clean up first####

rm(list=ls())

##### libraries for this script #####
library(cooccur)
library(tidyr)
library(ggplot2)
library(lmerTest)
library(lme4)
library(hrbrthemes)
library(RColorBrewer)
library(tidyverse)

#####set working directory & read in csv file#####
setwd("C:/Users/jansa/Desktop/2025-Mastomys-Final-Really")
setwd("C:/Users/jansa//Desktop/Ayo/2024_11_12_Post-Acacia-Final")
setwd("C:/Users/Jan_S/Desktop/2025-01-Mastomys-Final-Really")



M <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv", header = TRUE)


threshold <- 0.1*nrow(M) 

# Identify ManaMHCI_ columns
mana_cols <- grepl("^ManaMHC_", colnames(M))

# Calculate column sums for these columns
col_sums <- colSums(M[, mana_cols])

# Select columns with a sum greater than 53
keep_cols <- names(col_sums[col_sums > threshold])

# Create a new dataset with the filtered columns
M <- M[, c(setdiff(colnames(M), names(col_sums)), keep_cols)]







Mana <- subset(M, M$Village %in% c("Ebudin", "Aba gboro (Ile Ife)")) #Ayos 

Mana <- subset(M, M$Village %in% c("Ebudin"))

Mana <- subset(M, M$Village %in% c("Ekpoma", "Aba gboro (Ile Ife)", "Ebudin", "Owo")) #interesting nigeria

Mana <- subset(M, M$Country %in% c("Guinea")) # all guinea

Mana <- subset(M, M$Village %in% c("Damania", "Sokourala", "Ekpoma", "Ebudin", "Owo")) #positive nigeria and guinea


Mana <- subset(M, M$Village %in% c("Ekpoma", "Ebudin", "Owo")) #positive nigeria 
Mana <- subset(M, M$Village %in% c("Damania", "Sokourala", "Sonkonia")) # positive guinea

length(Mana$ID)
table(Mana$Village)


# Change 'neg' to 0 and 'pos' to 1 for the specified columns
Mana <- Mana %>%
  mutate(
    LASV_positive = ifelse(LASV_positive == "neg", 0, 1),
    IgG_positive = ifelse(IgG_positive == "neg", 0, 1)
  )
#### how many untrustworthy infection information pieces are in the data ####
#subset(bats, X229E_Correct.=="NO" | X2bBas_Correct.=="NO" | X2b_Correct.=="NO")
#### there are 13 #### -> remove them


#############        Co-occurence        ################

############################################################################
######### THE FOLLOWING IS FOR H.abae COOCCURANCE WITH SUPERTYPES#######
############################################################################

# the following file is a reduced version of my metadata file
# you need to have all columns you want to test for associations in binary (1=presence / 0=absence) format

#so we change the column in here but leave the data file unchanged
#lets do location first
#bats$AT <- ifelse(bats$Location=="Akpafu Todzi", 1, 0)
#bats$BUO <- ifelse(bats$Location=="Buoyem", 1, 0)
#bats$FO <- ifelse(bats$Location=="Forikrom", 1, 0)
#bats$KW <- ifelse(bats$Location=="Kwamang", 1, 0)
#bats$LT <- ifelse(bats$Location=="Likpe Todome", 1, 0)

#same with lineage/species
#bats$HcfruberB <- ifelse(bats$species=="H.cf.ruberB", 1, 0)
#bats$HcfruberC <- ifelse(bats$species=="H.cf.ruberC", 1, 0)
#bats$HcfruberD <- ifelse(bats$species=="H.cf.ruberD", 1, 0)
#bats$Habae <- ifelse(bats$species=="H.abae", 1, 0)

#lets also make a column for all totally uninfected as the inverse of the totally infected column

Mana$uninfected <- ifelse(Mana$number_alleles==1, 0, 1)

# co-occurrence analysis of common mhc alleles, STs and Lassavirus and IgG 
# we define common alleles as all mhc alleles that occur in at least 5 individuals or STs in 10

# first select all columns with allele and virus info we want to investigate for co-occurrence 
# either you select the columns by column number, or by column name (which is what I did)
#only consider STs common in more than 2 bats species (Kaesler et al.2017): ST6 was exluded as only present in C

Allelenames <- grep("^ManaMHC_", colnames(Mana), value = TRUE)

cooccur_Mana <- subset(Mana, select = c("LASV_positive", "IgG_positive", Allelenames))      # LASV, IgG, Supertypes




Supertypenames <- grep("^Supertype_", colnames(Mana), value = TRUE)

cooccur_Mana <- subset(Mana, select = c("LASV_positive", "IgG_positive", Supertypenames))   # Supertypes




# double-check there are no NAs in the matrix, because this will cause errors
# e.g. quick check all individuals are either 0 or 1 for allele ST1
# or if all individuals are either negative (0) or positive (1) for Coronavirus infection

cooccur_Mana$LASV_positive

head(cooccur_Mana) # have a quick look at the head of the table


# create a transposed matrix

cooccur_Mana_long <- t(as.matrix(cooccur_Mana))


# Filter out rows where the sum of the row values is 0
cooccur_Mana_long<- cooccur_Mana_long[rowSums(cooccur_Mana_long) != 0, ]

#Cooccur Mana 

# now apply the co-occurrence function
# just fyi when you look up the cooccurr package it shows co-occurrence of "species" and "sites",
# and these terms are also used within the function
# so here the "species" are alleles and pathogens and the "sites" are the individuals

ManaSTs <- cooccur(mat = cooccur_Mana_long, 
               type = "spp_site", 
               thresh = FALSE, 
               spp_names = TRUE)

  summary(ManaSTs) # short summary
plot(ManaSTs) # square plot with coloured positive / negative associations


# now we want to know which supertypes have a positive or negative association with lineages and infec stat

pair(mod = ManaSTs, "LASV_positive") 
pair(mod = ManaSTs, "IgG_positive") 

#Ni_ST_LASV_pair <-   pair(mod = ManaSTs, "LASV_positive") 
#Ni_ST_IgG_pair <- pair(mod = ManaSTs, "IgG_positive") 

Gu_ST_LASV_pair <-   pair(mod = ManaSTs, "LASV_positive") 
Gu_ST_IgG_pair <- pair(mod = ManaSTs, "IgG_positive") 



complete_ManaSTs<-prob.table(ManaSTs) # table of expected & observed cooccurrence and p-values for all comparisons
print(ManaSTs)
print(complete_ManaSTs)
obs.v.exp(ManaSTs) 

# effect.sizes(astro.alleles, standardized = FALSE, matrix = FALSE)       # raw effect sizes
effectsize_Mana  <- effect.sizes(ManaSTs, standardized = TRUE, matrix = FALSE)          # standardized effect sizes
complete_ManaSTs$effects <- effectsize_Mana$effects
complete_ManaSTs$obs.v.exp <- (complete_ManaSTs$exp_cooccur-complete_ManaSTs$obs_cooccur)*-1
# write.csv(complete_H.abae_STs, "complete_H.abae_STs.csv")  # you can write a csv file of standardized effect sizes

#heatmap.dat<-subset(heatmap.dat, sp2_name!="uninfected")


heatmap.dat<-subset(complete_ManaSTs, sp1_name=="LASV_positive" | sp1_name=="IgG_positive"  
                    #| sp1_name=="uninfected"
                    )
heatmap.dat<-subset(heatmap.dat, sp2_name!="IgG_positive")

#heatmap.dat<-subset(heatmap.dat, sp2_name!="uninfected")

heatmap.dat.p_gt <- heatmap.dat %>% 
  mutate(Asterisks_gt = ifelse(p_gt  <= 0.001, "***",
                            ifelse(p_gt <= 0.01, "**",
                                   ifelse(p_gt  <= 0.05, "*", NA))))
heatmap.dat.p_lt <- heatmap.dat %>% 
  mutate(Asterisks_lt = ifelse(p_lt  <= 0.001, "***",
                            ifelse(p_lt <= 0.01, "**",
                                   ifelse(p_lt  <= 0.05, "*", NA))))
#############
###### only one asteriks


heatmap.dat.p_gt <- heatmap.dat %>% 
  mutate(Asterisks_gt = ifelse(p_gt  <= 0.001, "*",
                               ifelse(p_gt <= 0.01, "*",
                                      ifelse(p_gt  <= 0.05, "*", NA))))
heatmap.dat.p_lt <- heatmap.dat %>% 
  mutate(Asterisks_lt = ifelse(p_lt  <= 0.001, "*",
                               ifelse(p_lt <= 0.01, "*",
                                      ifelse(p_lt  <= 0.05, "*", NA))))


############## HEATMAP SUPERTYPES 


ggplot(heatmap.dat, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  scale_fill_gradient2(high="red", mid="white", low="black", 
                       limits = c(-0.05, 0.05)
                       )+theme_minimal()+
 # scale_x_discrete(limits = c("Supertype_2","Supertype_3","Supertype_5","Supertype_10", "Supertype_13","Supertype_15",  "Supertype_19"))+
   scale_x_discrete(limits = c("Supertype_1","Supertype_2","Supertype_3",
                              "Supertype_4" , "Supertype_5","Supertype_6",
                              "Supertype_7", "Supertype_8", "Supertype_9",
                              "Supertype_10", "Supertype_11" ,"Supertype_12",
                              "Supertype_13","Supertype_14","Supertype_15", 
                              "Supertype_16" ,"Supertype_17", "Supertype_18",
                              "Supertype_19"))+
  theme(axis.text.x = element_text(angle = 90, size= 12))+
  theme(axis.text = element_text(size = 12))+
  ylab("Infection status")+xlab("Supertype")+
  geom_text(data = heatmap.dat.p_gt, aes(x = sp2_name, y = sp1_name, label = Asterisks_gt), size = 12)+
  geom_text(data = heatmap.dat.p_lt, aes(x = sp2_name, y = sp1_name, label = Asterisks_lt), size = 12)+
  coord_fixed(ratio = 1)


########################################### HEATMAP ALLELES
mana_columns <- grep("^ManaMHC_", colnames(Mana), value = TRUE)


ggplot(heatmap.dat, aes(sp2_name, sp1_name, fill=effects)) + 
  geom_tile()+
  scale_fill_gradient2(high="red", mid="white", low="darkgreen", 
                       limits = c(-0.05, 0.051)
  )+theme_minimal()+
  # scale_x_discrete(limits = c("Supertype_2","Supertype_3","Supertype_5","Supertype_10", "Supertype_13","Supertype_15",  "Supertype_19"))+
  scale_x_discrete(limits = c(mana_columns))+
  theme(axis.text.x = element_text(angle = 90))+
  ylab("Infection status")+xlab("Allele")+
  geom_text(data = heatmap.dat.p_gt, aes(x = sp2_name, y = sp1_name, label = Asterisks_gt), size = 4)+
  geom_text(data = heatmap.dat.p_lt, aes(x = sp2_name, y = sp1_name, label = Asterisks_lt), size = 4)+
  coord_fixed(ratio = 1)

#summary(glmer(CoV.229E ~ ST12 + (1|ID) , data=H.abae, family = binomial(link="logit"))) #confirms the results from the co-occurance 
#summary(glmer(uninfected ~ ST10 + (1|ID) , data=H.abae, family = binomial(link="logit"))) #same
