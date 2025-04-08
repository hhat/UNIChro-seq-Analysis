#!/usr/bin/env Rscript
# Bidirectional analysis functions

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
})

# Helper function for logit transformation
logit <- function(p) log(p/(1-p))

# Function to analyze a single SNP
analyze_snp <- function(snp_data) {
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
  model <- glmer(refalt ~ offset(logit(ALT_dna_prob)) + toALT_edit_bias + 
                 (1 + toALT_edit_bias | Donor),
                 family = binomial, data = long_DF)
  
  # Get fixed effects
  fixed_effects <- summary(model)$coefficients
  
  # Return results
  return(data.frame(
    SNP = unique(snp_data$SNP),
    effect = c("caQTL", "toALT_edit_bias"),
    Estimate = fixed_effects[,"Estimate"],
    Std_Error = fixed_effects[,"Std. Error"],
    p_value = fixed_effects[,"Pr(>|z|)"]
  ))
}

# Function to run analysis on full dataset
run_analysis <- function(data) {
  # Calculate ALT DNA probability
  result_DF <- data %>% mutate(ALT_dna_prob = round(ALT_count / (REF_count + ALT_count), digit=15)) 

  snps <- unique(result_DF$SNP)
  results <- data.frame()
  
  # Process each SNP
  for (snp in snps) {
    snp_data <- result_DF %>% filter(SNP == snp)
    snp_results <- analyze_snp(snp_data)
    results <- rbind(results, snp_results)
  }
  
  return(results)
}

# Function to read input data
read_input_data <- function(input_file) {
  read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

# Function to write results
write_results <- function(results, output_file) {
  write.table(results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

# Main function to run the entire pipeline
main <- function(input_file, output_file) {
  data <- read_input_data(input_file)
  results <- run_analysis(data)
  write_results(results, output_file)
  return(results)
}

# Check if script is being executed directly
script_name <- commandArgs()[grep("--file=", commandArgs(), fixed=TRUE)]
if (length(script_name) > 0) {  # Script is executed directly
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # Check for required arguments format
  if (length(args) < 4 || args[1] != "--input" || args[3] != "--output") {
    stop("Error: Invalid arguments.\nUsage: Rscript bidirectional_analysis.R --input input_file.txt --output results.txt")
  }
  
  # Get file paths
  input_file <- args[2]
  output_file <- args[4]
  
  # Verify input file exists
  if (!file.exists(input_file)) {
    stop("Error: Input file '", input_file, "' does not exist.")
  }
  
  # Run analysis
  main(input_file, output_file)
}
