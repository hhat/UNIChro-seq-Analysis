args <- commandArgs(trailingOnly = T)
IFILE1 <- as.character(args[1])  # REF allele UMI info
IFILE2 <- as.character(args[2])  # ALT allele UMI info
IFILE3 <- as.character(args[3])  # All UMIs from R1
IFILE4 <- as.character(args[4])  # UMIs from same barcode

edit_distance = 2
umi_length = 17

library(data.table)
library(magrittr)
library(dplyr)

OFILE1_prefix <- gsub(".tmp$", "", IFILE1)
OFILE2_prefix <- gsub(".tmp$", "", IFILE2)
CORRE_prefix <- gsub(".ref_umi_info.tmp$", "", IFILE1)

read_file <- function(file, snp_value) {
  if (file.exists(file)) {
    data <- fread(file)
    if (nrow(data) == 0) {
      return(NULL)
    }
    data <- as.data.frame(data)
    data$snp <- snp_value
    return(data)
  } else {
    return(NULL)
  }
}

d1_part1 <- read_file(IFILE1, "REF")
d1_part2 <- read_file(IFILE2, "ALT")

if (is.null(d1_part1) && is.null(d1_part2)) {
  stop("Both input files are empty or contain no rows.")
} else if (is.null(d1_part1)) {
  d1 <- d1_part2
} else if (is.null(d1_part2)) {
  d1 <- d1_part1
} else {
  d1 <- rbind(d1_part1, d1_part2)
}

colnames(d1) <- c("UMI", "R1_1st_base", "R2_1st_base", "CIGAR", "Size", "R2_Seq", "R2_strand", "snp")

TB <- table(d1$UMI)
umi_df <- data.frame(umi = names(TB), count = c(TB))

# Filter for valid UMIs (only ATGC)
is_valid_umi <- sapply(umi_df$umi, function(x) {
  x <- gsub("A", "", x)
  x <- gsub("T", "", x)
  x <- gsub("G", "", x)
  x <- gsub("C", "", x)
  nchar(x) == 0
})
umi_df <- umi_df[is_valid_umi, ]

# Filter for correct UMI length
is_valid_umi <- sapply(umi_df$umi, function(x) {
  nchar(x) == umi_length
})
umi_df <- umi_df[is_valid_umi, ]

# Sort by frequency (most common first)
umi_df <- umi_df[order(umi_df$count, decreasing = TRUE), ]

# Convert ATGC to binary representation
umi_df$umi_digit <- sapply(umi_df$umi, function(x) {
  x <- gsub("A", "1000", x)
  x <- gsub("T", "0100", x)
  x <- gsub("G", "0010", x)
  x <- gsub("C", "0001", x)
  return(x)
})

# Create UMI matrix for Hamming distance calculation
umi_matrix <- sapply(umi_df$umi_digit, function(x) {
  as.numeric(unlist(strsplit(split = "", x)))
})
colnames(umi_matrix) <- as.character(umi_df$umi)

# Calculate Hamming distances between UMIs
n_bases_matched <- t(umi_matrix) %*% umi_matrix
n_bases_mismatched <- umi_length - n_bases_matched

# Group similar UMIs using Hamming distance threshold
umi_df_mod <- umi_df
umi_corre <- data.frame(original_umi = character(), index_umi = character(), stringsAsFactors = FALSE)

if (nrow(umi_df_mod) == 1) {
  # Special case for single UMI
  single_umi <- umi_df_mod$umi[1]
  umi_corre <- data.frame(
    original_umi = single_umi,
    index_umi = single_umi,
    hamming_dist = 0,
    stringsAsFactors = FALSE
  )
} else {
  # Normal processing for multiple UMIs
  while (nrow(umi_df_mod) > 0) {
    # Get most frequent UMI as index
    index_umi <- umi_df_mod$umi[1]
    
    # Find UMIs within edit distance threshold
    diff_with_index <- n_bases_mismatched[index_umi, ]
    index_umi_cluster <- names(diff_with_index[diff_with_index <= edit_distance])
    
    # Only consider UMIs not yet assigned to clusters
    index_umi_cluster <- intersect(index_umi_cluster, umi_df_mod$umi)
    
    # Remove processed UMIs
    umi_df_mod <- umi_df_mod[!umi_df_mod$umi %in% index_umi_cluster, ]
    
    # Record UMI mapping
    umi_corre <- rbind(umi_corre, data.frame(
      original_umi = index_umi_cluster, 
      index_umi = rep(index_umi, length(index_umi_cluster)), 
      hamming_dist = as.numeric(n_bases_mismatched[index_umi, index_umi_cluster]), 
      stringsAsFactors = FALSE
    ))
  }
}

OFILE_corre <- paste0(CORRE_prefix, ".ED", edit_distance, ".UMI_correspondence.txt")
write.table(umi_corre, OFILE_corre, col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = "\t")

# Process all UMIs from sample fastq
fastq_umi <- fread(IFILE3, header = FALSE, data.table = FALSE)[, 1]

# Count UMI occurrences
fastq_umi_uniq <- as.matrix(table(fastq_umi))
fastq_umi_uniq <- data.frame(original_umi = rownames(fastq_umi_uniq), umi_count_fastq = fastq_umi_uniq)

# Apply UMI error correction
merged_df <- merge(fastq_umi_uniq, umi_corre, by = "original_umi", all.x = TRUE)
merged_df$index_umi <- ifelse(is.na(merged_df$index_umi), merged_df$original_umi, merged_df$index_umi)

# Summarize counts by index UMI
fastq_umi_uniq <- merged_df %>%
  group_by(index_umi) %>%
  summarise(umi_count_fastq = sum(umi_count_fastq)) %>% as.data.frame()

# Process UMIs from other samples with same barcode (for barcode hopping detection)
barcode_umi <- fread(IFILE4, header = FALSE, data.table = FALSE)[, 1]

# Count UMI occurrences
barcode_umi_uniq <- as.matrix(table(barcode_umi))
barcode_umi_uniq <- data.frame(original_umi = rownames(barcode_umi_uniq), umi_count_barcode = barcode_umi_uniq)

# Apply UMI error correction
merged_df <- merge(barcode_umi_uniq, umi_corre, by = "original_umi", all.x = TRUE)
merged_df$index_umi <- ifelse(is.na(merged_df$index_umi), merged_df$original_umi, merged_df$index_umi)

# Summarize counts by index UMI
barcode_umi_uniq <- merged_df %>%
  group_by(index_umi) %>%
  summarise(umi_count_barcode = sum(umi_count_barcode)) %>% as.data.frame()

# Process each allele (REF and ALT)
for (snp_value in c("REF", "ALT")) {
  if (snp_value == "REF") {
    d1_part <- d1_part1
    OFILE_prefix <- OFILE1_prefix
  } else if (snp_value == "ALT") {
    d1_part <- d1_part2
    OFILE_prefix <- OFILE2_prefix
  }

  if (!is.null(d1_part)) {
    colnames(d1_part) <- c("original_umi", "R1_1st_base", "R2_1st_base", "CIGAR", "Size", "R2_Seq", "R2_strand", "snp")
    d2 <- merge(d1_part, umi_corre, by = "original_umi")
    
    # Create position tags to identify PCR duplicates
    d2$original_ptag <- paste0(d2$original_umi, ":", d2$R1_1st_base, ":", d2$R2_1st_base,":",d2$Size)
    d2$ptag <- paste0(d2$index_umi, ":", d2$R1_1st_base, ":", d2$R2_1st_base,":",d2$Size)
    
    # Count reads per position tag
    d3 <- d2 %>%
      group_by(ptag) %>%
      mutate(ptag_count = n()) %>%
      ungroup() %>% as.data.frame()
    
    # Count total reads per UMI
    d4 <- d3 %>%
      group_by(index_umi) %>%
      mutate(count = n()) %>%
      ungroup() %>% as.data.frame()

    OFILE_data1 <- paste0(OFILE_prefix, ".ED", edit_distance, "_raw.txt")
    write.table(d4, OFILE_data1, col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = "\t")

    # Summarize UMI data
    d7 <- d4 %>%
      group_by(index_umi) %>%
      mutate(
        ptag_combined = paste(ptag, ptag_count, sep = "_count"),
        ptags = paste(names(sort(table(ptag_combined), decreasing = TRUE)), collapse = ","),
        Size_mode = as.numeric(names(sort(table(Size), decreasing = TRUE)[1]))
      ) %>%
      ungroup() %>%
      select(-hamming_dist, -original_umi, -original_ptag, -ptag, -ptag_count, -Size, -R1_1st_base, -R2_1st_base, -R2_Seq, -CIGAR, -ptag_combined) %>%
      distinct() %>%
      select(index_umi, everything()) %>%
      as.data.frame()
      
    # Add background UMI count from same sample
    d8 <- left_join(d7, fastq_umi_uniq)
    d8$samefastq_bg_woself_count <- d8$umi_count_fastq - d8$count 
    d8 <- d8 %>% select(-umi_count_fastq) %>% as.data.frame()
    
    # Add barcode hopping UMI count
    d9 <- left_join(d8, barcode_umi_uniq)
    d9[is.na(d9$umi_count_barcode), "umi_count_barcode"] <- 0

    # Write intermediate results
    OFILE_data1 <- paste0(OFILE_prefix, ".ED", edit_distance, "_index_umi_tmp.txt")
    write.table(d9, OFILE_data1, col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = "\t")
  }
}

# Finalize results and cross-reference REF/ALT counts
IFILE_1 <- paste0(OFILE1_prefix, ".ED", edit_distance, "_index_umi_tmp.txt")
IFILE_2 <- paste0(OFILE2_prefix, ".ED", edit_distance, "_index_umi_tmp.txt")
OFILE_1 <- paste0(OFILE1_prefix, ".ED", edit_distance, "_index_umi_final.txt") 
OFILE_2 <- paste0(OFILE2_prefix, ".ED", edit_distance, "_index_umi_final.txt") 

# Load processed data
d5REF <- if (file.exists(IFILE_1)) fread(IFILE_1, header = TRUE, data.table = FALSE) else data.frame()
d5ALT <- if (file.exists(IFILE_2)) fread(IFILE_2, header = TRUE, data.table = FALSE) else data.frame()

# Process REF data with ALT counts
if (nrow(d5REF) > 0) {
  if (nrow(d5ALT) > 0) {
    tmp <- distinct(d5ALT[, c("index_umi", "count")])
  } else {
    tmp <- data.frame(index_umi = character(), count = integer())
  }
  
  colnames(tmp)[2] <- "refaltchanged_count"
  REF_v2 <- left_join(d5REF, tmp, by = "index_umi")
  REF_v2$refaltchanged_count[is.na(REF_v2$refaltchanged_count)] <- 0
  
  REF_v3 <- REF_v2 %>% 
    select(index_umi, count, refaltchanged_count, samefastq_bg_woself_count, umi_count_barcode, everything())
} else {
  REF_v3 <- data.frame(
    index_umi = character(),
    count = integer(),
    refaltchanged_count = integer(),
    samefastq_bg_woself_count = integer(),
    umi_count_barcode = integer()
  )
}

# Process ALT data with REF counts
if (nrow(d5ALT) > 0) {
  if (nrow(d5REF) > 0) {
    tmp <- distinct(d5REF[, c("index_umi", "count")])
  } else {
    tmp <- data.frame(index_umi = character(), count = integer())
  }
  
  colnames(tmp)[2] <- "refaltchanged_count"
  ALT_v2 <- left_join(d5ALT, tmp, by = "index_umi")
  ALT_v2$refaltchanged_count[is.na(ALT_v2$refaltchanged_count)] <- 0
  
  ALT_v3 <- ALT_v2 %>% 
    select(index_umi, count, refaltchanged_count, samefastq_bg_woself_count, umi_count_barcode, everything())
} else {
  ALT_v3 <- data.frame(
    index_umi = character(),
    count = integer(),
    refaltchanged_count = integer(),
    samefastq_bg_woself_count = integer(),
    umi_count_barcode = integer()
  )
}

write.table(REF_v3, OFILE_1, col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = "\t")
write.table(ALT_v3, OFILE_2, col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE, sep = "\t")

# Clean up temporary files
file1 <- paste0(OFILE1_prefix, ".ED", edit_distance, "_index_umi_tmp.txt")
file2 <- paste0(OFILE2_prefix, ".ED", edit_distance, "_index_umi_tmp.txt")

if (file.exists(file1)) {
  file.remove(file1)
}

if (file.exists(file2)) {
  file.remove(file2)
}