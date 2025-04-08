#!/usr/bin/env Rscript

source("bidirectional_analysis.R")

run_analysis_perm <- function(data, seed = 12345) {
  set.seed(seed)
  
  result_DF <- data %>% mutate(ALT_dna_prob = round(ALT_count / (REF_count + ALT_count), digit=15)) 
  
  snps <- unique(result_DF$SNP)
  
  results <- data.frame()
  
  for(snp in snps) {
    snp_data <- result_DF %>% filter(SNP == snp)
    
    long_DF <- snp_data %>%
      select(SNP, Donor, REF_count, ALT_count, ALT_dna_prob, ref, alt, edit_direction) %>%
      pivot_longer(cols = c(ref, alt),
                 names_to = "refalt",
                 values_to = "count") %>%
      mutate(refalt = ifelse(refalt == "alt", 1, 0)) %>%
      uncount(count)
    
    long_DF <- long_DF %>%
      group_by(Donor, SNP) %>%
      mutate(
        edit_direction = sample(edit_direction),
        toALT_edit_bias = case_when(
          edit_direction == "REF_to_ALT" ~ 1,
          edit_direction == "NON_EDIT" ~ 0,
          edit_direction == "ALT_to_REF" ~ -1
        )
      ) %>%
      ungroup()
    
    # GLMM model
    model <- glmer(refalt ~ offset(logit(ALT_dna_prob)) + toALT_edit_bias + 
                  (1 + toALT_edit_bias | Donor),
                  family = binomial, data = long_DF,control = glmerControl(check.conv.singular = "ignore"))
    
    fixed_effects <- summary(model)$coefficients
    
    results <- rbind(results,
      data.frame(
        SNP = snp,
        effect_type = c("caQTL", "toALT_edit_bias"),
        Effect_size = fixed_effects[,"Estimate"],
        SE = fixed_effects[,"Std. Error"],
        p_value = fixed_effects[,"Pr(>|z|)"]
      )
    )
  }
  
  return(results)
}

run_permutation <- function(data, n_permutations = 100, seed = 123) {
  results_list <- list()
  base_seed <- seed
  
  # Setup progress display
  message("Starting permutation analysis with ", n_permutations, " iterations")
  
  for(i in 1:n_permutations) {
    
    iter_seed <- base_seed + i
    perm_results <- run_analysis_perm(data, seed = iter_seed)
    
    perm_results$permutation <- i
    results_list[[i]] <- perm_results
  }
  
  
  all_results <- bind_rows(results_list)
  return(all_results)
}

create_qq_data <- function(p_values) {
  n <- length(p_values)
  expected <- -log10((1:n) / (n+1))
  observed <- -log10(sort(p_values))
  
  return(data.frame(expected = expected, observed = observed))
}

plot_qq <- function(p_values, title = "QQ Plot", save_path = NULL) {
  qq_data <- create_qq_data(p_values)
  
  max_val <- max(c(qq_data$expected, qq_data$observed)) * 1.05
  
  plot <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = title,
      x = "Expected -log10(p)",
      y = "Observed -log10(p)"
    ) + 
    theme_minimal() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      strip.text = element_text(size = 20),
      plot.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      axis.line = element_line(color = "black"),
      panel.grid.minor = element_blank(),  
      panel.grid.major = element_blank()
    ) +
    xlim(0, max_val) + 
    ylim(0, max_val)
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot, width = 7, height = 6)
  }
  
  return(plot)
}