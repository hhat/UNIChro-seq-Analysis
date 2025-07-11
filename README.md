# UNIChro-seq-Analysis

This repository provides analysis scripts used in the UNIChro-seq paper (in preparation).

## Folder Structure

- **`Unichro-seq-pipeline/`**  
  - **`ATAC-seq/`**  
    - **`script/`**: Contains scripts for preprocessing and analyzing ATAC-seq data processing.
    - `workflow.sh`: Workflow script for ATAC-seq analysis.
  - **`DNA-seq/`**  
    - **`script/`**: Contains scripts for preprocessing and analyzing DNA-seq data processing.
    - `workflow.sh`: Workflow script for DNA-seq analysis.
      
- **`crisper_analysis/`**  
  - **`demo/`**: Contains a demo of detecting a caQTL effect using unidirectional editing.
    - `README.md`: Documentation for the demo.
    - `uni-directional_analysis.R`: R script for .
    - `input_file_REF_to_ALT.txt`: Input file example.
    - `results_REF_to_ALT.txt`: Analysis result example.
    - `input_file_ALT_to_REF.txt`: Input file example.
    - `results_ALT_to_REF.txt`: Analysis result example.
  - **`CRISPResso2/`**: Contains scripts for analysis using CRISPResso2.
  
- **`others/`**
  - **`siLEF1_RNA-seq/`**: Contains scripts for siLEF1 RNA-seq data analysis.
  - **`siLEF1_ATAC-seq/`**: Contains scripts for siLEF1 ATAC-seq data analysis.

## uni-directional Editing Demo

The [`crisper_analysis/uni-directional_edit_demo/`](./crisper_analysis/uni-directional_edit_demo/) directory includes a demonstration of uni-directional editing analysis. For detailed information, please refer to the `README.md` within the respective directory.

