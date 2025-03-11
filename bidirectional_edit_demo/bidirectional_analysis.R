#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
})

logit <- function(p) log(p/(1-p))

# Command-line argument parsing
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript glm_analysis.R --input input_file.txt --output output_file.txt")
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
  select(Donor, SNP, ref, alt, REF_count, ALT_count, edit_direction) %>% 
  mutate(ALT_dna_prob = ALT_count / (REF_count + ALT_count))

# Unique SNPs
snps <- unique(result_DF$SNP)
results <- data.frame()

# Process each SNP
for (snp in snps) {
  snp_data <- result_DF %>% filter(SNP == snp)

  # Prepare data in long format
      long_DF <- snp_data %>%
      select(SNP, Donor, REF_count, ALT_count, ALT_dna_prob, ref, alt, edit_direction) %>%
      pivot_longer(cols = c(ref, alt),
                   names_to = "refalt",
                   values_to = "count") %>%
      mutate(refalt = ifelse(refalt == "alt", 1, 0)) %>%
      uncount(count) %>% 
      mutate(
        toALT_edit_bias = case_when(
          edit_direction == "REF_to_ALT" ~ 1,
          edit_direction == "NON_EDIT" ~ 0,
          edit_direction == "ALT_to_REF" ~ -1
        )
      )

  # GLMM with random slope
  model <- glmer(refalt ~ offset(logit(ALT_dna_prob)) + toALT_edit_bias + (1 + toALT_edit_bias | Donor),
                          family = binomial, data = long_DF)
   fixed_effects <- summary(model)$coefficients
  # Store results
  results <- rbind(results,
    data.frame(
      SNP = snp,
      effect = c("caQTL", "toALT_edit_bias"),
      Estimate = fixed_effects[,"Estimate"],
      Std_Error = fixed_effects[,"Std. Error"],
      p_value = fixed_effects[,"Pr(>|z|)"]
    )
  )

}

# Save results
write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("[INFO] GLMM analysis completed. Results saved to:", output_file, "\n")
