#!/usr/bin/env Rscript

# -------------------------------------------------------------------------
# Script Name: calc_delta_significance.R
# Description: 
#   Calculates the Delta statistic between Unconditional and Conditional 
#   GWAS results to distinguish growth-specific signals from plumage-confounded ones.
#   
# Input: 
#   1. Unconditional GWAS summary statistics (GEMMA output)
#   2. Conditional GWAS summary statistics (GEMMA output with covariate)
#
# Output:
#   A tab-separated file containing merged statistics, Delta values, 
#   and effect categories for each SNP.
# -------------------------------------------------------------------------

library(tidyverse) # Load core tidyverse packages

# =========================================================================
# 1. Load Summary Statistics
# =========================================================================
# Model 1: Unconditional (y = W + x + u + e)
# Model 2: Conditional (y = W + c + x + u + e), where c is plumage color

message("Loading unconditional GWAS results...")
# read_table automatically handles whitespace delimiters common in GEMMA outputs.
# col_types = cols() suppresses the column specification message for cleaner logs.
model1 <- read_table("output/gwas_unconditional.assoc.txt", col_types = cols())

message("Loading conditional GWAS results...")
model2 <- read_table("output/gwas_conditional.assoc.txt", col_types = cols())

# =========================================================================
# 2. Data Merging
# =========================================================================
# Merge datasets by SNP ID ('rs'). 
# inner_join ensures only SNPs present in both analyses are retained.
# 'suffix' appends distinct labels to columns with identical names (e.g., p_lrt).
merged <- inner_join(model1, model2, by = "rs", suffix = c("_uncond", "_cond"))

# =========================================================================
# 3. Calculate Delta and Categorize Effects
# =========================================================================
# Formula: Delta = -log10(P_cond) - [-log10(P_uncond)]
# Interpretation:
#   Delta > 0: Significance increased after conditioning (Signal revealed/strengthened)
#   Delta < 0: Significance decreased after conditioning (Signal driven by confounder)

# Define epsilon threshold for "negligible change"
# This threshold effectively filters out noise around zero change.
epsilon <- 0.1 

final_results <- merged %>%
  mutate(
    # Transform P-values to -log10 scale
    logP_uncond = -log10(p_lrt_uncond),
    logP_cond   = -log10(p_lrt_cond),
    
    # Calculate the Delta statistic
    Delta = logP_cond - logP_uncond,
    
    # Categorize SNPs based on Delta direction and magnitude
    Effect_Category = case_when(
      Delta > epsilon  ~ "Increased Significance (Growth Signal)",
      Delta < -epsilon ~ "Reduced Significance (Confounded by Plumage)",
      TRUE             ~ "Stable Signal" # Captures remaining cases where abs(Delta) <= epsilon
    )
  )

# =========================================================================
# 4. Save Results
# =========================================================================
output_file <- "output/paired_gwas_delta_results.txt"
message(paste("Saving paired GWAS Delta results to:", output_file))

# write_tsv saves as tab-delimited text, standard for bioinformatics data
write_tsv(final_results, output_file)

message("Processing complete. Summary of effect categories:")
print(table(final_results$Effect_Category))

# For Reproducibility: Print session info
sessionInfo()
