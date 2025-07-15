# UNIChro-seq-Analysis

[![DOI](https://zenodo.org/badge/800069394.svg)](https://doi.org/10.5281/zenodo.15892171)
This repository provides analysis scripts used in the UNIChro-seq paper (in preparation).

## Folder Structure

* `Unichro-seq-pipeline/`
   * `ATAC-seq/`
      * `script/`: Contains scripts for preprocessing and analyzing ATAC-seq data.
      * `workflow.sh`: Workflow script for ATAC-seq analysis.
   * `DNA-seq/`
      * `script/`: Contains scripts for preprocessing and analyzing DNA-seq data.
      * `workflow.sh`: Workflow script for DNA-seq analysis.
   * `design_primer/`
      * `script/`: Contains scripts for primer design and evaluation.
      * `workflow.sh`: Workflow script for primer design.
   * `crispr_analysis/`
      * `uni-directional_edit_demo/`: Contains a demo for detecting caQTL effects using uni-directional editing in homozygous samples (REF/REF or ALT/ALT).      
      * `README.md`: Documentation for the demo.
      * `uni-directional_analysis.R`: R script for main analysis.
      * `uni-directional_editing_demo.ipynb`: Jupyter notebook for the demo.
      * `input_file_REF_to_ALT.txt`: Input file example.
      * `results_REF_to_ALT.txt`: Analysis result example.
      * `input_file_ALT_to_REF.txt`: Input file example2.
      * `results_ALT_to_REF.txt`: Analysis result example2.
      * `images/`: Contains concept images for the demo.
   * `CRISPResso2/`: Contains scripts for analysis using CRISPResso2.
   * `others/`
   * `siLEF1_RNA-seq/`: Contains scripts for siLEF1 RNA-seq data analysis.
   * `siLEF1_ATAC-seq/`: Contains scripts for siLEF1 ATAC-seq data analysis.

## Uni-directional Editing Demo

The [`crispr_analysis/uni-directional_edit_demo/`](crispr_analysis/uni-directional_edit_demo/) directory includes a demonstration of uni-directional editing analysis for homozygous samples (REF/REF or ALT/ALT). For detailed information, please refer to the `README.md` within the respective directory.

## License

This project is licensed under the terms of the LICENSE file included in this repository.
