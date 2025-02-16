# UNIChro-seq-Analysis

## Analysis of Bi-directional Editing
## Identification of caQTL Effects and Validation of Edit Bias

## Project Purpose
This project aims to identify **chromatin accessibility QTLs (caQTLs)** and verify that observed allele-specific accessibility is not due to technical artifacts (referred to as “edit bias”). We combine data from **probe-based DNA-seq** and **probe-based ATAC-seq**, both targeting the same genomic sites, to evaluate how genetic variants (e.g., SNPs) influence chromatin accessibility.

## Data Used
- **Probe-based DNA-seq Count Data**  
  Counts of reference and alternate alleles obtained via targeted DNA sequencing using probes. This serves as a baseline for the true allele frequency.
- **Probe-based ATAC-seq Count Data**  
  Counts of reference and alternate alleles from ATAC-seq (which measures chromatin accessibility), collected with the same probe-based approach for the same genomic regions as DNA-seq. By comparing both datasets, we can assess whether accessibility differs significantly by allele.

All count data are located in the `data/` directory (see the **Demo Dataset Structure** below for details).

## Analysis Steps

1. **Data Preprocessing**  
   - Merge DNA-seq and ATAC-seq count data. For each SNP (or genomic locus), compile the reference and alternate allele reads.  
   - Filter out low-coverage or poor-quality entries as needed.

2. **GLM Modeling**  
   - Use a **generalized linear model (GLM)** for each SNP to test for allele-specific accessibility.  
   - For instance, compare the ratio of reference vs. alternate reads in ATAC-seq against the expected 50:50 ratio from DNA-seq using logistic regression.

3. **Statistical Significance**  
   - Extract effect sizes (e.g., odds ratios, log odds) and p-values from the GLM, then apply multiple-testing correction.  
   - Sites passing the significance threshold are flagged as potential caQTLs.

4. **Meta-Analysis** (if required)  
   - If multiple samples or conditions exist, aggregate effect sizes across datasets using meta-analysis (e.g., a random-effects model).  
   - This approach accounts for between-sample variability and produces a combined effect estimate.

5. **Edit Bias Validation**  
   - For each significant caQTL, verify that DNA-seq does not exhibit a similar allele imbalance.  
   - If DNA-seq and ATAC-seq both show the same allele bias, it may indicate a technical artifact (edit bias). Those sites can be excluded or re-evaluated.

## Output Results
- **List of Detected caQTLs**  
  A table (e.g., `caQTL_effects.tsv`) with effect sizes, p-values, and significance indicators for each SNP.
- **Forest Plot**  
  A visualization (e.g., `forest_plot.png`) showing effect sizes and confidence intervals across multiple samples/conditions, and the combined meta-analysis result.
- **Edit Bias Analysis**  
  A table (e.g., `edit_bias_analysis.tsv`) comparing DNA-seq and ATAC-seq allele frequencies to confirm no inherent bias.  
  - If DNA-seq reads are near a 50:50 split but ATAC-seq shows a strong imbalance, it suggests a true allele-specific accessibility signal.

All such results are typically placed in the `results/` directory.

## Demo Dataset Structure
A possible repository structure for this demo is:

