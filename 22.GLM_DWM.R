library(lubridate)
library(lme4)
library(ggplot2)
library(dplyr)
setwd("C:/Users/Jan_S/Desktop/2025-01-Mastomys-Final-Really")


load("22.GLM.RData")

mana_mhc <- read.csv("Nigeria_Guinea_Wide_ST_noGaps.csv")

str(mana_mhc)
mana_mhc <- mana_mhc %>%
  mutate(year_capture = year(ymd(Date_Capture)))

plyr::ddply(mana_mhc, c("Country"), summarise, n=length(Country))
#Country   n
#1  Guinea 315
#2 Nigeria 350
plyr::ddply(mana_mhc, c("Country", "Village", "IgG_positive", "LASV_positive"), summarise, n=length(Country))

#### unique alleles Nigeria 
nigeria_mhc<-subset(mana_mhc, Country=="Nigeria")
guinea_mhc<-subset(mana_mhc, Country=="Guinea")

# Sum across rows for each column
col_sums_nig <- colSums(nigeria_mhc[,17:409])
col_sums_gui <- colSums(guinea_mhc[,17:409])

# Logical vectors of presence/absence for each allele
present_nig <- col_sums_nig >= 1
present_gui <- col_sums_gui >= 1

# Alleles unique to Nigeria = present in Nigeria but absent in Guinea
unique_nig <- present_nig & !present_gui
sum(unique_nig)   # number of unique alleles in Nigeria

# Alleles unique to Guinea = present in Guinea but absent in Nigeria
unique_gui <- present_gui & !present_nig
sum(unique_gui)   # number of unique alleles in Guinea

# If you want actual allele names (column names):
#names(unique_nig[unique_nig])
#names(unique_gui[unique_gui])

# Alleles shared in Guinea & in Nigeria
shared <- present_gui & present_nig
sum(shared)   # number of unique alleles in Guinea



ggplot(mana_mhc, aes(number_alleles))+geom_histogram(bins = 52)
plyr::ddply(mana_mhc, c("Country"), summarise, mean(number_alleles), sd(number_alleles))
t.test(subset(mana_mhc, Country=="Nigeria")$number_alleles, subset(mana_mhc, Country=="Guinea")$number_alleles)
plyr::ddply(mana_mhc, c("Country"), summarise, mean(number_ST), sd(number_ST))
t.test(subset(mana_mhc, Country=="Nigeria")$number_ST, subset(mana_mhc, Country=="Guinea")$number_ST)

mana_mhc$LASV_positive  <- as.factor(mana_mhc$LASV_positive)

mana_mhc$IgG_positive   <- as.factor(mana_mhc$IgG_positive)

mana_mhc$IgG_positive[mana_mhc$IgG_positive == "."] <- NA


mana_mhc$year_capture   <- as.factor(mana_mhc$year_capture )
mana_mhc$Village   <- as.factor(mana_mhc$Village )



subset_NA <- mana_mhc[is.na(mana_mhc$Sex) | is.na(mana_mhc$ELW) | is.na(mana_mhc$IgG_positive), ]

subset_NA[, c("ID", "Village", "Sex", "ELW", "IgG_positive" )]

#mana_mhc_subset <- subset(mana_mhc, Village!="Tambaya" & Village!="Aba gboro (Ile Ife)" & Village!="Ifon" & Village!="Okeluse" & Village!="Okhuesan")

# Drop rows with NA in any predictor
model_df <- mana_mhc %>%
  filter(!is.na(ELW), !is.na(Sex), !is.na(IgG_positive))



model_df <- model_df %>%
  mutate(Sex = case_when(
    Sex %in% c("M", "m") ~ "M",
    Sex %in% c("F", "f") ~ "F",
    TRUE ~ as.character(Sex)
  )) %>%
  mutate(Sex = factor(Sex))  # re-factor after correction






######

library(glmmTMB)
lm2_LassaPCR<-(glmmTMB(LASV_positive ~ Supertype_5+Country + #here the ManaMHCI.006 name is needed
                         #number_ST +
                         ELW + 
                         #Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm2_LassaPCR)
summary(lm2_LassaPCR)


lm2_LassaPCR<-(glmmTMB(LASV_positive ~ Supertype_18*Country + #here the ManaMHCI.006 name is needed
                         #number_ST +
                         ELW + 
                         #Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm2_LassaPCR)
summary(lm2_LassaPCR)

lm2_LassaPCR<-(glmmTMB(IgG_positive ~ Supertype_5+Country + #here the ManaMHCI.006 name is needed
                         #number_ST +
                         ELW + 
                         Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm2_LassaPCR)
summary(lm2_LassaPCR)

lm2_LassaPCR<-(glmmTMB(IgG_positive ~ Supertype_15*Country + #here the ManaMHCI.006 name is needed
                         number_ST +
                         ELW + 
                         Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm2_LassaPCR)
summary(lm2_LassaPCR)

################
alleles_LASV <- c("ManaMHC_017", "ManaMHC_048", "ManaMHC_104")
alleles_IgG <- c("ManaMHC_009", "ManaMHC_011", "ManaMHC_012", "ManaMHC_021", 
                 "ManaMHC_022", "ManaMHC_025", "ManaMHC_027", "ManaMHC_033", 
                 "ManaMHC_050", "ManaMHC_069", "ManaMHC_104", "ManaMHC_107",
                 "ManaMHC_197")

###### LASV loop 


results_LASV <- list()

for (al in alleles_LASV) {
  formula_lasv <- as.formula(
    paste("LASV_positive ~", al, "* Country + number_alleles + 
          ELW + 
          Sex +
          (1|year_capture)")
  )
  model <- glmer(formula_lasv, data = model_df, family = binomial())
  results_LASV[[al]] <- summary(model)
}

results_LASV


# IGg loop

results_igg <- list()

for (al in alleles_IgG) {
  formula_igg <- as.formula(
    paste("IgG_positive ~", al, "* Country + number_alleles + 
          ELW + 
          Sex + 
          (1|year_capture)")
  )
  model <- glmer(formula_igg, data = model_df, family = binomial())
  results_igg[[al]] <- summary(model)
}

results_igg
############## FDR correction #####

# LASV
p_values_lasv_glmer <- data.frame(
  Term = character(),
  Allele = character(),
  Estimate = numeric(), # Added Estimate column
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (al in names(results_LASV)) {
  summary_al <- results_LASV[[al]]$coefficients
  allele_name <- gsub("ManaMHC_", "", al)
  
  for (i in 1:nrow(summary_al)) {
    term_name <- rownames(summary_al)[i]
    estimate_val <- summary_al[i, "Estimate"] # Extract the Estimate
    p_value <- summary_al[i, "Pr(>|z|)"]
    
    p_values_lasv_glmer <- rbind(p_values_lasv_glmer, data.frame(
      Term = term_name,
      Allele = allele_name,
      Estimate = estimate_val, # Include the Estimate
      P_Value = p_value
    ))
  }
}

# Perform FDR correction on the p-values
p_values_lasv_glmer_fdr <- p_values_lasv_glmer %>%
  group_by(Allele) %>%
  mutate(FDR_P_Value = p.adjust(P_Value, method = "fdr")) %>%
  ungroup()

print(p_values_lasv_glmer_fdr)
print(p_values_lasv_glmer_fdr, n = nrow(p_values_lasv_glmer_fdr))

significant_lasv_fdr <- p_values_lasv_glmer_fdr %>%
  filter(FDR_P_Value < 0.05)

print(significant_lasv_fdr)
print(significant_lasv_fdr, n = nrow(significant_lasv_fdr))


####### model + plot Mana 017

#### model looking into the negative association between acute LASV infection and 006 ####
library(glmmTMB)
lm1_LassaPCR<-(glmmTMB(LASV_positive ~ ManaMHC_017*Country + 
                         number_alleles +
                         ELW + 
                         Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm1_LassaPCR)
summary(lm1_LassaPCR)

lm2_LassaPCR<-(glmmTMB(LASV_positive ~ ManaMHC_048+Country + 
                         number_alleles +
                         ELW + 
                         Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm2_LassaPCR)
summary(lm2_LassaPCR)

lm3_LassaPCR<-(glmmTMB(LASV_positive ~ ManaMHC_104 + Country +
                         number_alleles +
                         ELW + 
                         Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm3_LassaPCR)
summary(lm3_LassaPCR)

lm4_LassaPCR<-(glmmTMB(LASV_positive ~ Supertype_5 * Country +
                         #number_ST +
                         ELW + 
                         Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm4_LassaPCR)
summary(lm4_LassaPCR)

lm5_LassaPCR<-(glmmTMB(LASV_positive ~ Supertype_18 * Country +
                         number_ST +
                         ELW + 
                         Sex + 
                         (1|year_capture), 
                       data=model_df, family=binomial()))
MuMIn::dredge(lm5_LassaPCR)
summary(lm5_LassaPCR)

#######

( allsets <- emmeans(lm3_LassaPCR, ".") )
contrast(allsets, "eff")


# your contrast results (example subset for ManaMHC_009)
contrast_df <- data.frame(
  contrast = c("ManaMHC_1040", 
               "ManaMHC_1041"),
  estimate = c(-0.352, 0.352),
  SE = c(0.191, 0.191)
)

# split into two columns
contrast_df$ManaMHC_104 <- ifelse(grepl("ManaMHC_1040", contrast_df$contrast), 0, 1)

# 95% CI
contrast_df$lower <- contrast_df$estimate - 1.96 * contrast_df$SE
contrast_df$upper <- contrast_df$estimate + 1.96 * contrast_df$SE

#contrast_df[1,4]<-0

# plot
LASV104_plot<-ggplot(contrast_df, aes(x=as.factor(ManaMHC_104), y=estimate)) +
  geom_point(position=position_dodge(width=0.4), size=4) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2,
                position=position_dodge(width=0.4)) +
  geom_line(position=position_dodge(width=0.4)) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  labs(y="LASV+ (estimate ± 95% CI)", 
       x="Mana*104") +ylim(-1.5,1.5)+
  theme_minimal(base_size = 14)+theme(axis.title = element_text(face="bold"), legend.position = "none")

#######

( allsets <- emmeans(lm1_LassaPCR, ".") )
contrast(allsets, "eff")


# your contrast results 
contrast_df <- data.frame(
  contrast = c("ManaMHC_0170 Guinea", "ManaMHC_0171 Guinea",
               "ManaMHC_0170 Nigeria", "ManaMHC_0171 Nigeria"),
  estimate = c(0.414, -0.552, -0.635, 0.773),
  SE = c(0.229, 0.466, 0.295, 0.261)
)

# split into two columns
contrast_df$ManaMHC_017 <- ifelse(grepl("ManaMHC_0170", contrast_df$contrast), 0, 1)
contrast_df$Country <- ifelse(grepl("Guinea", contrast_df$contrast), "Guinea", "Nigeria")

# 95% CI
contrast_df$lower <- contrast_df$estimate - 1.96 * contrast_df$SE
contrast_df$upper <- contrast_df$estimate + 1.96 * contrast_df$SE


# plot
LASV17_plot<-ggplot(contrast_df, aes(x=as.factor(ManaMHC_017), y=estimate, color=Country, group=Country)) +
  geom_point(position=position_dodge(width=0.4), size=4) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2,
                position=position_dodge(width=0.4)) +
  geom_line(position=position_dodge(width=0.4)) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  labs(y="LASV+ (estimate ± 95% CI)", 
       x="Mana*017",
       color="Country") +ylim(-1.5,1.5)+
  scale_colour_manual(values=c("#D29D01", "#00CC7D"))+
  theme_minimal(base_size = 14)+theme(axis.title = element_text(face="bold"), legend.position = "none")




###### IgG FDR

p_values_igg_glmer <- data.frame(
  Term = character(),
  Allele = character(),
  Estimate = numeric(), # Added Estimate column
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (al in names(results_igg)) {
  summary_al <- results_igg[[al]]$coefficients
  allele_name <- gsub("ManaMHC_", "", al)
  
  for (i in 1:nrow(summary_al)) {
    term_name <- rownames(summary_al)[i]
    estimate_val <- summary_al[i, "Estimate"] # Extract the Estimate
    p_value <- summary_al[i, "Pr(>|z|)"]
    
    p_values_igg_glmer <- rbind(p_values_igg_glmer, data.frame(
      Term = term_name,
      Allele = allele_name,
      Estimate = estimate_val, # Include the Estimate
      P_Value = p_value
    ))
  }
}

# Perform FDR correction on the p-values
p_values_igg_glmer_fdr <- p_values_igg_glmer %>%
  group_by(Allele) %>%
  mutate(FDR_P_Value = p.adjust(P_Value, method = "fdr")) %>%
  ungroup()

print(p_values_igg_glmer_fdr)
print(p_values_igg_glmer_fdr, n = nrow(p_values_igg_glmer_fdr))

significant_igg_fdr <- p_values_igg_glmer_fdr %>%
  filter(FDR_P_Value < 0.05)

print(significant_igg_fdr)
print(significant_igg_fdr, n = nrow(significant_igg_fdr))


############### IgG Mana 009 model and plot

IgG009<-(glmer(IgG_positive ~ ManaMHC_009 * Country + #here the ManaMHCI.006 name is needed
                number_alleles +
                ELW + 
                Sex + 
                (1|year_capture), 
              data=model_df, family=binomial()))
MuMIn::dredge(IgG009)
summary(IgG009)

( allsets <- emmeans(IgG009, ".") )
contrast(allsets, "eff")


# your contrast results (example subset for ManaMHC_009)
contrast_df <- data.frame(
  contrast = c("ManaMHC_0090", 
               "ManaMHC_0091"),
  estimate = c(-0.523, 0.523),
  SE = c(0.161, 0.161)
)

# split into two columns
contrast_df$ManaMHC_009 <- ifelse(grepl("ManaMHC_0690", contrast_df$contrast), 0, 1)

# 95% CI
contrast_df$lower <- contrast_df$estimate - 1.96 * contrast_df$SE
contrast_df$upper <- contrast_df$estimate + 1.96 * contrast_df$SE

contrast_df[1,4]<-0

# plot
IgG009_plot<-ggplot(contrast_df, aes(x=as.factor(ManaMHC_009), y=estimate)) +
  geom_point(position=position_dodge(width=0.4), size=4) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2,
                position=position_dodge(width=0.4)) +
  geom_line(position=position_dodge(width=0.4)) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  labs(y="IgG+ (estimate ± 95% CI)", 
       x="Mana*009") +ylim(-1.5,1.5)+
  theme_minimal(base_size = 14)+theme(axis.title = element_text(face="bold"), legend.position = "none")




############### IgG Mana 69 model and plot

IgG69<-(glmer(IgG_positive ~ ManaMHC_069 * Country + #here the ManaMHCI.006 name is needed
                       number_alleles +
                       ELW + 
                       Sex + 
                       (1|year_capture), 
                     data=model_df, family=binomial()))
MuMIn::dredge(IgG69)
summary(IgG69)

( allsets <- emmeans(IgG69, ".") )
contrast(allsets, "eff")


# your contrast results (example subset for ManaMHC_069 × Country)
contrast_df <- data.frame(
  contrast = c("ManaMHC_0690 Guinea", "ManaMHC_0691 Guinea",
               "ManaMHC_0690 Nigeria", "ManaMHC_0691 Nigeria"),
  estimate = c(0.491, -0.328, -0.264, 0.101),
  SE = c(0.218, 0.286, 0.212, 0.322)
)

# split into two columns
contrast_df$ManaMHC_069 <- ifelse(grepl("ManaMHC_0690", contrast_df$contrast), 0, 1)
contrast_df$Country <- ifelse(grepl("Guinea", contrast_df$contrast), "Guinea", "Nigeria")

# 95% CI
contrast_df$lower <- contrast_df$estimate - 1.96 * contrast_df$SE
contrast_df$upper <- contrast_df$estimate + 1.96 * contrast_df$SE


# plot
IgG69_plot<-ggplot(contrast_df, aes(x=as.factor(ManaMHC_069), y=estimate, color=Country, group=Country)) +
  geom_point(position=position_dodge(width=0.4), size=4) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2,
                position=position_dodge(width=0.4)) +
  geom_line(position=position_dodge(width=0.4)) +
  geom_hline(yintercept = 0, linetype="dashed", color="black") +
  labs(y="IgG+ (estimate ± 95% CI)", 
       x="Mana*069",
       color="Country") + ylim(-1.5,1.5)+
  scale_colour_manual(values=c("#D29D01", "#00CC7D"))+
  theme_minimal(base_size = 14)+theme(axis.title = element_text(face="bold"), legend.position = "none")


ggpubr::ggarrange(LASV104_plot, LASV17_plot, IgG009_plot, IgG69_plot, nrow=2, ncol=2, labels=c("A", "B", "C", "D"))


##### supertypes IgG

IgGST5<-(glmer(IgG_positive ~ Supertype_5 * Country + #here the ManaMHCI.006 name is needed
                 number_ST +
                 ELW + 
                 Sex + 
                 (1|year_capture), 
               data=model_df, family=binomial()))
MuMIn::dredge(IgGST5)
summary(IgGST5)

IgGST15<-(glmer(IgG_positive ~ Supertype_15 * Country + #here the ManaMHCI.006 name is needed
                 number_ST +
                 ELW + 
                 Sex + 
                 (1|year_capture), 
               data=model_df, family=binomial()))
MuMIn::dredge(IgGST15)
summary(IgGST15)

( allsets <- emmeans(IgGST15, ".") )
contrast(allsets, "eff")

###### EWL
############### individual Tests EWL
library(glmmTMB)
lm_ELW09<-(glmmTMB((ELW) ~ ManaMHC_009*IgG_positive*Country+
                       number_alleles +
                       Sex + 
                       (1|year_capture), 
                     data=model_df, family=gaussian()))
MuMIn::dredge(lm_ELW09)
shapiro.test(resid(lm_ELW09))
summary(lm_ELW09)
ggplot(model_df, aes(IgG_positive, ELW))+geom_boxplot()
ggplot(model_df, aes(IgG_positive, ELW, fill=as.factor(ManaMHC_009)))+geom_boxplot()+facet_grid(~Country)


ggplot(model_df, aes(IgG_positive, ELW, fill=as.factor(ManaMHC_069)))+geom_boxplot()+facet_grid(~Country)+
  theme_minimal()+theme(axis.title = element_text(size=14, face="bold"), legend.position = c(.1,.88))+
  labs(y="eye lens weight (mg)", 
       x="IgG detection",
       fill="Mana*069")+scale_fill_manual(values=c("grey90", "grey55"))

( allsets <- emmeans(lm_ELW69, ".") )
contrast(allsets, "eff")




