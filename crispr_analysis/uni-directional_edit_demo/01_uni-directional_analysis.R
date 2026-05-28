#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

logit <- function(p) log(p/(1-p))

analyze_snp <- function(snp_data) {
  long_DF <- snp_data %>%
    select(Donor, SNP ,NON_EDITED_DNA_count, EDITED_DNA_count, non_edited_atac_count, edited_atac_count) %>%
    mutate(EDITED_DNA_prob = round(EDITED_DNA_count / (NON_EDITED_DNA_count + EDITED_DNA_count), digit=15)) %>%
    pivot_longer(cols = c(non_edited_atac_count, edited_atac_count),
                 names_to = "edit_type",
                 values_to = "count") %>%
    mutate(edit_type = ifelse(edit_type == "edited_atac_count", 1, 0)) %>%
    uncount(count)
  
  # GLMM with random slope
  model <- glmer(edit_type ~ offset(logit(EDITED_DNA_prob)) + (1 | Donor),
                 family = binomial, data = long_DF)
  
  fixed_effects <- summary(model)$coefficients
  
  return(data.frame(
    SNP = unique(long_DF$SNP),
    Effect_size = fixed_effects[1,"Estimate"],
    SE = fixed_effects[1,"Std. Error"],
    p_value = fixed_effects[1,"Pr(>|z|)"]
  ))
}

run_analysis <- function(data) {

  snps <- unique(data$SNP)
  results <- data.frame()
  
  for (snp in snps) {
    snp_data <- data %>% filter(SNP == snp)
    snp_results <- analyze_snp(snp_data)
    results <- rbind(results, snp_results)
  }
  
  return(results)
}

read_input_data <- function(input_file) {
  read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

write_results <- function(results, output_file) {
  write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

main <- function(input_file, output_file) {
  data <- read_input_data(input_file)
  results <- run_analysis(data)
  write_results(results, output_file)
  return(results)
}

# Check whether script is executed directly
script_name <- commandArgs()[grep("--file=", commandArgs(), fixed=TRUE)]
if (length(script_name) > 0) {  # Script is executed directly
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check for required arguments format
  if (length(args) < 4 || args[1] != "--input" || args[3] != "--output") {
    stop("Error: Invalid arguments.\nUsage: Rscript bidirectional_analysis.R --input input_file --output output_file")
  }
  
  input_file <- args[2]
  output_file <- args[4]
  
  main(input_file, output_file)
}
