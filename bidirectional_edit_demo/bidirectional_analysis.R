#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# Command-line argument parsing
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript glm_analysis.R --input input_file.tsv --output output_file.tsv")
}

input_file <- output_file <- NULL

for (i in seq_along(args)) {
  if (args[i] == "--input") {
    input_file <- args[i + 1]
  } else if (args[i] == "--output") {
    output_file <- args[i + 1]
  }
}

cat("[INFO] Input file  :", input_file, "\n")
cat("[INFO] Output file :", output_file, "\n")

# Read input data (selecting only required columns)
result_DF <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(sample_id, SNP, ref, alt, REF_count, ALT_count, edit_direction)

# Unique SNPs
snps <- unique(result_DF$SNP)
results_glm <- data.frame()

# Process each SNP
for (snp in snps) {
  snp_data <- result_DF %>% filter(SNP == snp)

  # Prepare data in long format
  long_DF <- snp_data %>%
    select(SNP, REF_count, ALT_count, ref, alt, edit_direction) %>%
    pivot_longer(cols = c(ref, alt), names_to = "refalt", values_to = "count") %>%
    mutate(refalt = ifelse(refalt == "alt", 1, 0)) %>%
    uncount(count) %>%
    mutate(
      toALT_edit_linear = case_when(
        edit_direction == "REF_to_ALT" ~ 1,
        edit_direction == "NON_EDIT" ~ 0,
        edit_direction == "ALT_to_REF" ~ -1
      ),
      toALT_edit_nonlinear = ifelse(edit_direction == "REF_to_ALT", 1, 0)
    )

  # GLM fitting
  model <- glm(refalt ~ offset(logit(ALT_count / (REF_count + ALT_count))) + 
               toALT_edit_linear + toALT_edit_nonlinear,
               family = binomial, data = long_DF)
  summary_model <- summary(model)

  # Store results
  results_glm <- rbind(results_glm, data.frame(
    SNP = snp,
    effect = "Intercept",
    Estimate = coef(summary_model)["(Intercept)", "Estimate"],
    Std_Error = coef(summary_model)["(Intercept)", "Std. Error"],
    p_value = coef(summary_model)["(Intercept)", "Pr(>|z|)"]
  ))

  results_glm <- rbind(results_glm, data.frame(
    SNP = snp,
    effect = "toALT_edit_linear",
    Estimate = coef(summary_model)["toALT_edit_linear", "Estimate"],
    Std_Error = coef(summary_model)["toALT_edit_linear", "Std. Error"],
    p_value = coef(summary_model)["toALT_edit_linear", "Pr(>|z|)"]
  ))

  results_glm <- rbind(results_glm, data.frame(
    SNP = snp,
    effect = "toALT_edit_nonlinear",
    Estimate = coef(summary_model)["toALT_edit_nonlinear", "Estimate"],
    Std_Error = coef(summary_model)["toALT_edit_nonlinear", "Std. Error"],
    p_value = coef(summary_model)["toALT_edit_nonlinear", "Pr(>|z|)"]
  ))
}

# Save results
write.table(results_glm, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO] GLM analysis completed. Results saved to:", output_file, "\n")
