# README: GLM with an Offset Using DNA-seq Ratio

This repository provides an R script that uses **DNA-seq allele counts** (with columns `REF, ALT`) as an offset in a logistic regression, allowing for the detection of **allele-specific accessibility biases** in **ATAC-seq** data (with columns `ref, alt`). 

---

## Script Overview

- **Script Name:** `glm_modeling_offset.R`
- **Goal:**
  1. Account for DNA-seq allele ratios (which might not be 50:50) by using them as an **offset**.
  2. Perform logistic regression on ATAC-seq (`ref, alt`) counts to estimate **caQTL effects** (allelic imbalance) relative to the DNA-seq baseline.

- **Column Naming Convention:**
  - **DNA-seq:** `SNP, REF, ALT` (uppercase)
  - **ATAC-seq:** `SNP, ref, alt` (lowercase)
  - Both contain a common column, `SNP`, used for joining the data frames.

- **GLM Model (Conceptual):**  
  \[
    \mathrm{cbind}( alt_{\mathrm{ATAC}},\ ref_{\mathrm{ATAC}} ) 
    \sim 1\ +\ \mathrm{offset}\Bigl(\logit\Bigl(\frac{ALT_{\mathrm{DNA}}}{REF_{\mathrm{DNA}} + ALT_{\mathrm{DNA}}}\Bigr)\Bigr)
  \]  
  - Here, the **intercept** captures how much the ATAC-seq ratio deviates from the DNA-seq ratio.

---

## Input Data Example

### DNA-seq (`dna_counts.tsv`)

| SNP   | REF  | ALT  |
|-------|------|------|
| rs001 | 100  | 130  |
| rs002 |  80  |  90  |
| ...   | ...  | ...  |

- `REF`: Reference allele read count  
- `ALT`: Alternate allele read count  

### ATAC-seq (`atac_counts.tsv`)

| SNP   | ref  | alt  |
|-------|------|------|
| rs001 |  50  |  60  |
| rs002 |  40  |  50  |
| ...   | ...  | ...  |

- `ref`: Reference allele read count  
- `alt`: Alternate allele read count  

---

## How to Run the Script

1. **Install Dependencies**

   - **R** (version 3.6 or later recommended)
   - **CRAN Packages**: `dplyr`, `tidyr`
     ```r
     install.packages(c("dplyr", "tidyr"))
     ```

2. **Execute the Script**
   ```bash
   Rscript bidirectional_edit_analysis.R \
       --dna data/dna_counts.tsv \
       --atac data/atac_counts.tsv \
       --output results/glm_results.tsv
