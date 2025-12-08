library(ggplot2)
library(dplyr)

# Assuming your glmer model summaries for LASV are in the list 'results_LASV'

load("23.Forestplots.RData")

forest_data_lasv <- data.frame()

alleles_to_plot <- c("ManaMHC_017", "ManaMHC_048", "ManaMHC_104")

interaction_terms <- paste0("ManaMHC_", sprintf("%03d", as.numeric(gsub("ManaMHC_", "", alleles_to_plot))), ":CountryNigeria")
terms_to_plot <- c(alleles_to_plot, interaction_terms)

for (term in terms_to_plot) {
  allele_name <- gsub("ManaMHC_", "", term)
  allele_part <- gsub(":CountryNigeria", "", term)
  summary_model <- NULL
  if (grepl(":", term)) {
    if (!is.null(results_LASV[[allele_part]])) {
      summary_model <- results_LASV[[allele_part]]$coefficients
    }
  } else {
    if (!is.null(results_LASV[[term]])) {
      summary_model <- results_LASV[[term]]$coefficients
    }
  }
  
  if (!is.null(summary_model)) {
    term_row <- grep(paste0("^", gsub(":", "\\:", term), "$"), rownames(summary_model))
    
    if (length(term_row) > 0) {
      effect_size <- summary_model[term_row, "Estimate"]
      std_error <- summary_model[term_row, "Std. Error"]
      lower_ci <- effect_size - 1.96 * std_error # Approximate 95% CI
      upper_ci <- effect_size + 1.96 * std_error
      
      fdr_p_value_df <- p_values_lasv_glmer_fdr %>%
        filter(Allele == gsub("ManaMHC_", "", allele_part), Term == term) %>%
        pull(FDR_P_Value)
      
      significance <- ifelse(!is.na(fdr_p_value_df) && fdr_p_value_df < 0.05, "*", "")
      
      forest_data_lasv <- rbind(forest_data_lasv, data.frame(
        Term = term,
        Effect = effect_size,
        LowerCI = lower_ci,
        UpperCI = upper_ci,
        Significance = significance
      ))
    } else {
      cat("Warning: Term", term, "not found in model summary\n")
    }
  } else {
    cat("Warning: No results found for term", term, "\n")
  }
}

# Create the forest plot
ggplot(forest_data_lasv, aes(y = Term, x = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.2) +
  geom_point(size = 3) +
  geom_text(aes(label = Significance, y = Term, x = Effect), vjust = -1, hjust = 0.5, size = 10) +
  labs(title = "Effect of Alleles and Allele:Country Interactions on LASV Positivity (FDR Corrected)",
       x = "Log Odds Ratio (Effect Size)",
       y = "") +
  theme_minimal()

###########################    IGG #########

library(ggplot2)
library(dplyr)

# Assuming your glmer model summaries for IgG are in the list 'results_igg'
# and your FDR-adjusted p-values for IgG are in 'p_values_igg_glmer_fdr'

forest_data_igg <- data.frame()

# Define the alleles you want to plot
alleles_igg_to_plot <- unique(gsub("ManaMHC_", "", p_values_igg_glmer_fdr$Allele))

interaction_terms_igg <- paste0("ManaMHC_", sprintf("%03d", as.numeric(alleles_igg_to_plot)), ":CountryNigeria")
terms_igg_to_plot <- c(paste0("ManaMHC_", sprintf("%03d", as.numeric(alleles_igg_to_plot))), interaction_terms_igg)

for (term in terms_igg_to_plot) {
  allele_name <- gsub("ManaMHC_", "", term)
  allele_part <- gsub(":CountryNigeria", "", term)
  summary_model <- NULL
  if (grepl(":", term)) {
    if (!is.null(results_igg[[allele_part]])) {
      summary_model <- results_igg[[allele_part]]$coefficients
    }
  } else {
    if (!is.null(results_igg[[term]])) {
      summary_model <- results_igg[[term]]$coefficients
    }
  }
  
  if (!is.null(summary_model)) {
    term_row <- grep(paste0("^", gsub(":", "\\:", term), "$"), rownames(summary_model))
    
    if (length(term_row) > 0) {
      effect_size <- summary_model[term_row, "Estimate"]
      std_error <- summary_model[term_row, "Std. Error"]
      lower_ci <- effect_size - 1.96 * std_error # Approximate 95% CI
      upper_ci <- effect_size + 1.96 * std_error
      
      fdr_p_value_df <- p_values_igg_glmer_fdr %>%
        filter(Allele == allele_name, Term == term) %>%
        pull(FDR_P_Value)
      
      significance <- ifelse(!is.na(fdr_p_value_df) && fdr_p_value_df < 0.05, "*", "")
      error_color <- ifelse(effect_size > 0, "red", "darkgreen")
      
      forest_data_igg <- rbind(forest_data_igg, data.frame(
        Term = term,
        Effect = effect_size,
        LowerCI = lower_ci,
        UpperCI = upper_ci,
        Significance = significance,
        ErrorColor = error_color
      ))
    } else {
      cat("Warning: Term", term, "not found in model summary\n")
    }
  } else {
    cat("Warning: No results found for term", term, "\n")
  }
}

# Create the forest plot for IgG with colored error bars
ggplot(forest_data_igg, aes(y = Term, x = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI, color = ErrorColor), height = 0.2) +
  scale_color_identity() + # Use the color defined in the data frame
  geom_point(size = 3) +
  geom_text(aes(label = Significance, y = Term, x = Effect), vjust = -.01, hjust = 0.5, size = 10) +
  theme(axis.text.y = element_text(size = 20)) +
  labs(title = "Effect of Alleles and Allele:Country Interactions on IgG Presence (FDR Corrected)",
       x = "Log Odds Ratio (Effect Size)",
       y = "") +
  theme_minimal()



############# IGG Only Alleles, no interaction forest plot:##############

results_IGg <- list()

for (al in alleles_IgG) {
  formula_igg <- as.formula(
    paste("IgG_positive ~", al, "* Country + number_alleles + ELW + Sex + (1|year_capture)")
  )
  model <- glmer(formula_igg, data = model_df, family = binomial())
  results_IGg[[al]] <- model  # store full model, not summary
}

# Assuming p_values_igg_glmer_fdr already exists (use same pipeline as for LASV but on results_IGg)

forest_df_main_alleles_igg <- data.frame(
  Term = character(),
  Allele = character(),
  Estimate = numeric(),
  Conf.low = numeric(),
  Conf.high = numeric(),
  P_Value = numeric(),
  FDR_P_Value = numeric(),
  stringsAsFactors = FALSE
)

for (al in names(results_IGg)) {
  model <- results_IGg[[al]]
  model_summary <- summary(model)$coefficients
  allele_name <- gsub("ManaMHC_", "", al)
  
  if (al %in% rownames(model_summary)) {
    # Estimate, p-value
    est <- model_summary[al, "Estimate"]
    p_val <- model_summary[al, "Pr(>|z|)"]
    
    # Wald confidence intervals
    conf_ints <- confint(model, parm = al, method = "Wald")
    ci_low <- conf_ints[1]
    ci_high <- conf_ints[2]
    
    # FDR-corrected p-value
    fdr_val <- p_values_igg_glmer_fdr %>%
      filter(Term == al, Allele == allele_name) %>%
      pull(FDR_P_Value)
    
    # Append to dataframe
    forest_df_main_alleles_igg <- rbind(forest_df_main_alleles_igg, data.frame(
      Term = al,
      Allele = allele_name,
      Estimate = est,
      Conf.low = ci_low,
      Conf.high = ci_high,
      P_Value = p_val,
      FDR_P_Value = fdr_val
    ))
  }
}



####################
# Prepare data: add significance markers and color coding based on effect direction
forest_df_main_alleles_igg <- forest_df_main_alleles_igg %>%
  mutate(
    Significance = ifelse(FDR_P_Value < 0.05, "*", ""),
    ErrorColor = ifelse(Estimate > 0, "red", "darkgreen")
  )

# Plot
ggplot(forest_df_main_alleles_igg, aes(y = Term, x = Estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = Conf.low, xmax = Conf.high, color = ErrorColor), height = 0.2) +
  scale_color_identity() +  # Use ErrorColor values directly
  geom_point(size = 3) +
  geom_text(aes(label = Significance), vjust = -.01, hjust = 0.5, size = 10) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 20)) +
  labs(
    title = "Effect of Alleles on IgG Presence (FDR Corrected)",
    x = "Log Odds Ratio (Effect Size)",
    y = ""
  )




#####


# Clean allele names and sort by allele number
forest_df_main_alleles_igg <- forest_df_main_alleles_igg %>%
  mutate(
    AlleleNum = as.numeric(gsub("ManaMHC_", "", Term)),   # extract numeric part
    TermClean = paste0("*", sprintf("%03d", AlleleNum)),  # format as *050, *104, etc.
    Significance = ifelse(FDR_P_Value < 0.05, "*", ""),
    ErrorColor = ifelse(Estimate > 0, "red", "darkgreen")
  ) %>%
  arrange(AlleleNum)

# Plot
ggplot(forest_df_main_alleles_igg, aes(y = reorder(TermClean, AlleleNum), x = Estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = Conf.low, xmax = Conf.high, color = ErrorColor), height = 0.2) +
  scale_color_identity() +  # Use ErrorColor values directly
  geom_point(size = 3) +
  geom_text(aes(label = Significance), vjust = -.01, hjust = 0.5, size = 10) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 20)) +
  labs(
    title = "Effect of Alleles on IgG Presence (FDR Corrected)",
    x = "Log Odds Ratio (Effect Size)",
    y = ""
  )




################ Country Interaction of 069

# Define the terms to extract from the model summary
# The effect for Guinea is the main effect of the allele
term_guinea <- "ManaMHC_069"
# The interaction term is the difference in effect for Nigeria
term_nigeria_interaction <- "ManaMHC_069:CountryNigeria"

# Initialize an empty dataframe for forest plot data
forest_data_069_interaction <- data.frame()

# Get the full model summary for allele 069
summary_model <- results_igg[["ManaMHC_069"]]$coefficients

# Check if the summary exists and contains the necessary terms
if (!is.null(summary_model) && term_guinea %in% rownames(summary_model) && term_nigeria_interaction %in% rownames(summary_model)) {
  
  # 1. Extract and calculate data for Guinea
  effect_guinea <- summary_model[term_guinea, "Estimate"]
  std_error_guinea <- summary_model[term_guinea, "Std. Error"]
  lower_ci_guinea <- effect_guinea - 1.96 * std_error_guinea
  upper_ci_guinea <- effect_guinea + 1.96 * std_error_guinea
  
  # Get the FDR p-value for the main effect of the allele
  fdr_p_value_guinea <- p_values_igg_glmer_fdr %>%
    filter(Allele == "069", Term == term_guinea) %>%
    pull(FDR_P_Value)
  
  significance_guinea <- ifelse(!is.na(fdr_p_value_guinea) && fdr_p_value_guinea < 0.05, "*", "")
  
  forest_data_069_interaction <- rbind(forest_data_069_interaction, data.frame(
    Term = "ManaMHC_069 in Guinea",
    Effect = effect_guinea,
    LowerCI = lower_ci_guinea,
    UpperCI = upper_ci_guinea,
    Significance = significance_guinea,
    ErrorColor = ifelse(effect_guinea > 0, "red", "darkgreen")
  ))
  
  # 2. Extract and calculate data for Nigeria (main effect + interaction effect)
  effect_nigeria_interaction <- summary_model[term_nigeria_interaction, "Estimate"]
  
  # Calculate the combined effect and standard error for Nigeria
  # This requires the variance-covariance matrix to correctly calculate the standard error of the sum
  v_cov_matrix <- as.matrix(vcov(results_igg[["ManaMHC_069"]]))
  se_nigeria <- sqrt(v_cov_matrix[term_guinea, term_guinea] + v_cov_matrix[term_nigeria_interaction, term_nigeria_interaction] + 2 * v_cov_matrix[term_guinea, term_nigeria_interaction])
  
  effect_nigeria <- effect_guinea + effect_nigeria_interaction
  lower_ci_nigeria <- effect_nigeria - 1.96 * se_nigeria
  upper_ci_nigeria <- effect_nigeria + 1.96 * se_nigeria
  
  # Get the FDR p-value for the interaction term
  fdr_p_value_nigeria <- p_values_igg_glmer_fdr %>%
    filter(Allele == "069", Term == term_nigeria_interaction) %>%
    pull(FDR_P_Value)
  
  significance_nigeria <- ifelse(!is.na(fdr_p_value_nigeria) && fdr_p_value_nigeria < 0.05, "*", "")
  
  forest_data_069_interaction <- rbind(forest_data_069_interaction, data.frame(
    Term = "ManaMHC_069 in Nigeria",
    Effect = effect_nigeria,
    LowerCI = lower_ci_nigeria,
    UpperCI = upper_ci_nigeria,
    Significance = significance_nigeria,
    ErrorColor = ifelse(effect_nigeria > 0, "red", "darkgreen")
  ))
}

# The ggplot code will now work because forest_data_069_interaction is populated
forest_data_069_interaction$Country <- ifelse(grepl("Nigeria", forest_data_069_interaction$Term), "Nigeria", "Guinea")

ggplot(forest_data_069_interaction, aes(y = Term, x = Effect, color = Country)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI, color = Country), height = 0.2) +
  geom_point(size = 3) +
  geom_text(aes(label = Significance, y = Term, x = Effect), vjust = -0.5, size = 10) +
  scale_color_manual(values = c("darkgreen", "red")) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 15)) +
  labs(
    title = "Effect of ManaMHC_069 on IgG Presence by Country",
    x = "Log Odds Ratio (Effect Size)",
    y = ""
  )


################### ELW 
library(dplyr)
library(tibble)
library(ggplot2)
library(lme4) # Ensure lme4 is loaded for the lmer models

# Assuming 'results_elw' is populated with lmer model objects
# And 'elw_p_values_all_effects' is populated with FDR corrected p-values for ELW models

forest_df_elw_plot <- data.frame(
  Term = character(),
  Allele = character(),
  Label = character(), # Custom label for the Y-axis
  Estimate = numeric(),
  Conf.low = numeric(),
  Conf.high = numeric(),
  P_Value = numeric(),
  FDR_P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each allele's model results
for (al_col_name in names(results_elw)) {
  model <- results_elw[[al_col_name]]
  model_summary <- summary(model)$coefficients
  allele_number <- gsub("ManaMHC_", "", al_col_name) # Extract the allele number (e.g., "009")
  
  # --- Include all IgG_positivepos main effects for each allele ---
  if ("IgG_positivepos" %in% rownames(model_summary)) {
    term_name <- "IgG_positivepos"
    est <- model_summary[term_name, "Estimate"]
    std_err <- model_summary[term_name, "Std. Error"] # Get standard error
    p_val <- model_summary[term_name, "Pr(>|t|)"]
    
    # --- Calculate Wald confidence intervals manually to avoid 'long vectors' error ---
    # Using a z-score (1.96 for 95% CI) as a common approximation for large samples/models
    ci_low <- est - 1.96 * std_err
    ci_high <- est + 1.96 * std_err
    
    # Get the FDR-corrected p-value from our pre-calculated table
    fdr_val <- elw_p_values_all_effects %>%
      filter(Term == term_name, Allele == allele_number) %>%
      pull(FDR_P_Value)
    
    # Add to our plotting dataframe
    forest_df_elw_plot <- rbind(forest_df_elw_plot, data.frame(
      Term = term_name,
      Allele = allele_number,
      Label = paste0("IgG+ (Allele ", allele_number, ")"), # Clear label for plot
      Estimate = est,
      Conf.low = ci_low,
      Conf.high = ci_high,
      P_Value = p_val,
      FDR_P_Value = fdr_val
    ))
  }
  
  # --- Include only the SIGNIFICANT three-way interaction: ManaMHC_069:IgG_positivepos:CountryNigeria ---
  # This interaction term only exists for allele 069
  if (allele_number == "069") { 
    three_way_term <- paste0("ManaMHC_", allele_number, ":IgG_positivepos:CountryNigeria")
    if (three_way_term %in% rownames(model_summary)) {
      # Check if this specific interaction term is significant based on its FDR_P_Value
      fdr_val_interaction <- elw_p_values_all_effects %>%
        filter(Term == three_way_term, Allele == allele_number) %>%
        pull(FDR_P_Value)
      
      # Only add to the plot if it's significant (FDR_P_Value < 0.05)
      if (!is.null(fdr_val_interaction) && length(fdr_val_interaction) > 0 && fdr_val_interaction < 0.05) {
        est <- model_summary[three_way_term, "Estimate"]
        std_err <- model_summary[three_way_term, "Std. Error"] # Get standard error
        p_val <- model_summary[three_way_term, "Pr(>|t|)"]
        
        # --- Calculate Wald confidence intervals manually ---
        ci_low <- est - 1.96 * std_err
        ci_high <- est + 1.96 * std_err
        
        forest_df_elw_plot <- rbind(forest_df_elw_plot, data.frame(
          Term = three_way_term,
          Allele = allele_number,
          Label = paste0("Allele ", allele_number, ":IgG+:CountryInteraction"), # Clear label for plot
          Estimate = est,
          Conf.low = ci_low,
          Conf.high = ci_high,
          P_Value = p_val,
          FDR_P_Value = fdr_val_interaction
        ))
      }
    }
  }
}

# --- Prepare data for plotting: add significance markers and color coding ---
forest_df_elw_plot <- forest_df_elw_plot %>%
  mutate(
    Significance = ifelse(FDR_P_Value < 0.05, "*", ""), # Add '*' for significant terms
    ErrorColor = ifelse(Estimate > 0, "red", "darkgreen") # Color based on positive/negative effect
  ) %>%
  # Order the terms for the Y-axis: IgG+ main effects first, then interaction, ordered by allele number
  arrange(ifelse(grepl("IgG_positivepos$", Term), 1, 2), # Puts 'IgG_positivepos' terms before interactions
          as.numeric(Allele)) # Then orders numerically by allele

# --- Create the Forest Plot using ggplot2 ---
ggplot(forest_df_elw_plot, aes(y = reorder(Label, -as.numeric(factor(Label, levels = Label))), x = Estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + # Vertical line at zero effect
  geom_errorbarh(aes(xmin = Conf.low, xmax = Conf.high, color = ErrorColor), height = 0.2) + # Horizontal error bars
  scale_color_identity() + # Use the defined colors directly
  geom_point(size = 3) + # Point for the estimate
  geom_text(aes(label = Significance), vjust = -0.01, hjust = 0.5, size = 10) + # Significance asterisks
  theme_minimal() + # Minimal theme for clean look
  theme(axis.text.y = element_text(size = 12)) + # Adjust Y-axis label size for readability
  labs(
    title = "Effect of IgG Status and Interaction on ELW (FDR Corrected)", # Plot title
    x = "Coefficient (Effect Size on ELW)", # X-axis label (ELW is a continuous outcome)
    y = "" # Y-axis label (labels come from 'Label' aesthetic)
  )

