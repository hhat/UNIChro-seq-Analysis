#!/usr/bin/env Rscript

# Load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set filtering parameters
count_cutoff <- as.numeric(10)

# Initialize results dataframe
res <- data.frame()

# Process batches
for (batch in c("20231117")) {
  # Set base directory for this batch
  base_dir <- paste0("/home/imgkono/wd/img/crispr_qtl/UNIChro_seq2/bowtie/", batch, "/")
  
  # Find all bowtie directories
  bowtie_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  
  # Initialize files vector
  files <- c()
  
  # Find all filtered files for the batch
  for (dir in bowtie_dirs) {
    files <- c(files, list.files(path = dir, 
                pattern = paste0("ED2_index_umi_final.txt_size300_count10_halffilter.txt"),
                recursive = TRUE, 
                full.names = TRUE))
  }
  
  # Process each file and collect row counts
  for (ite in 1:length(files)) {
    file <- files[ite]
    DFumi <- read.table(file, header = TRUE, sep = "\t")
    res <- rbind(res, data.frame(filename = file, nrow = nrow(DFumi)))
  }
}

# Initialize columns for metadata
res$sequence <- rep(NA, nrow(res))
res$ID <- rep(NA, nrow(res))
res$barcode <- rep(NA, nrow(res))
res$target <- rep(NA, nrow(res))
res$snp <- rep(NA, nrow(res))
res$size_cutoff <- rep(NA, nrow(res))
res$count_cutoff <- rep(NA, nrow(res))

# Extract metadata from filenames
for (ite in 1:nrow(res)) {
  filename <- res[ite, 1]
  # Normalize path separators
  normalized_filename <- gsub("//+", "/", filename)
  
  # Use regex to extract components from filename
  regex <- ".*/([0-9]{8}(?:_[0-9]{8})?)/([^/]+)\\.([^/]+)/([^/]+)/\\4\\.(ref|alt)_.*_size([0-9]+)_count([0-9]+)_.*"
  matches <- regmatches(normalized_filename, regexec(regex, normalized_filename))
  
  # Store extracted components
  result <- matches[[1]][-1]  
  
  res$sequence[ite] <- result[1]        
  res$ID[ite] <- result[2]      
  res$barcode[ite] <- result[3]     
  res$target[ite] <- result[4]      
  res$snp[ite] <- result[5]    
  res$size_cutoff[ite] <- result[6]        
  res$count_cutoff[ite] <- result[7]                    
}

# Remove filename column and reshape data
fin <- res[, -1]
data <- fin %>%
  pivot_wider(names_from = snp, values_from = nrow) 

# Rename columns and sort
colnames(data)[4] <- "SNP"
data <- data[order(data$SNP), ]
data$sequence <- as.numeric(data$sequence)

# Write output
write.csv(data, "processed_umi_counts.csv", row.names = FALSE)