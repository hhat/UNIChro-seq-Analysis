# UNIChro-seq Bidirectional Analysis Demo

This demo-script identifies caQTLs and editing biases using a generalized linear model that combines UMI-counted ATAC-seq reads and DNA-seq ALT/(REF+ALT) ratios, both analyzed with allele-specific probes.

## Usage

### Running the Script
```bash
Rscript bidirectional_analysis.R --input data/input_file.txt --output results/results.txt
```

### Command-line Arguments
| Argument | Description |
|----------|-------------|
| `--input` | Path to the input text file |
| `--output` | Path to the output text file |

## Expected Input Format (`input_file.txt`)
The script expects an tab-sepqrated input file with the following required columns:

| Donor | edit_direction | SNP | ref | alt | REF_count | ALT_count |
|--------|---------------|-----|-----|-----|-----------|-----------| 
| R11_01 | ALT_to_REF | rs2248137 | 1930 | 963 | 111093 | 51338 |
| R11_01 | NON_EDIT | rs2248137 | 1303 | 1423 | 76226 | 57463 |
| R11_01 | REF_to_ALT | rs2248137 | 1512 | 2016 | 68108 | 74336 |
| R11_02 | ALT_to_REF | rs2248137 | 917 | 527 | 127240 | 55054 |
| R11_02 | NON_EDIT | rs2248137 | 1016 | 908 | 99598 | 83489 |
| R11_02 | REF_to_ALT | rs2248137 | 900 | 1188 | 65335 | 71255 |
| R11_03 | ALT_to_REF | rs2248137 | 1002 | 707 | 90922 | 46997 |
| R11_03 | NON_EDIT | rs2248137 | 763 | 721 | 99369 | 75986 |
| R11_03 | REF_to_ALT | rs2248137 | 683 | 1329 | 85982 | 76401 |
| R11_04 | ALT_to_REF | rs2248137 | 980 | 502 | 50745 | 25753 |
| R11_04 | NON_EDIT | rs2248137 | 800 | 733 | 37124 | 34223 |
| R11_04 | REF_to_ALT | rs2248137 | 920 | 959 | 38982 | 52319 |
| R11_05 | ALT_to_REF | rs2248137 | 866 | 476 | 41834 | 22113 |
| R11_05 | NON_EDIT | rs2248137 | 684 | 678 | 30891 | 29560 |
| R11_05 | REF_to_ALT | rs2248137 | 693 | 927 | 33109 | 42735 |
| R11_06 | ALT_to_REF | rs2248137 | 1373 | 713 | 37531 | 19821 |
| R11_06 | NON_EDIT | rs2248137 | 677 | 1115 | 33244 | 31417 |
| R11_06 | REF_to_ALT | rs2248137 | 802 | 1024 | 30610 | 39797 |
| R11_07 | ALT_to_REF | rs2248137 | 924 | 545 | 44439 | 22701 |
| R11_07 | NON_EDIT | rs2248137 | 663 | 670 | 39037 | 35905 |
| R11_07 | REF_to_ALT | rs2248137 | 862 | 1339 | 29703 | 39208 |
| R11_08 | ALT_to_REF | rs2248137 | 810 | 636 | 37079 | 20072 |
| R11_08 | NON_EDIT | rs2248137 | 453 | 623 | 31634 | 29393 |
| R11_08 | REF_to_ALT | rs2248137 | 790 | 1628 | 24162 | 33027 |
| R11_09 | ALT_to_REF | rs2248137 | 590 | 496 | 38228 | 21528 |
| R11_09 | NON_EDIT | rs2248137 | 279 | 324 | 32131 | 29908 |
| R11_09 | REF_to_ALT | rs2248137 | 418 | 747 | 28405 | 36058 |

### Column Descriptions
* `Donor`: Unique identifier for the Donor
* `edit_direction`: Editing direction (`REF_to_ALT`, `ALT_to_REF`, or `NON_EDIT`)
* `SNP`: SNP ID (chromosome and position)
* `ref`: Reference allele read count from ATAC-seq
* `alt`: Alternative allele read count from ATAC-seq
* `REF_count`: Reference allele read count from DNA-seq
* `ALT_count`: Alternative allele read count from DNA-seq

## Output Format (`results.txt`)
The script generates an output file containing the results of the GLMM analysis in the following format:

| SNP | effect | Estimate | Std_Error | p_value |
|-----|---------|-----------|------------|----------|
| rs2248137 | caQTL | 0.2006 | 0.0528 | 0.000144 |
| rs2248137 | toALT_edit_bias | 0.0169 | 0.0316 | 0.592 |

### Output Column Descriptions
* `SNP`: The SNP ID being analyzed
* `effect`: The type of effect
  * `caQTL`: caQTL effect (positive values indicate higher chromatin accessibility for ALT allele)
  * `toALT_edit_bias`: Edit bias effect (positive values indicate bias towards ALT allele)
* `Estimate`: The estimated coefficient value
* `Std_Error`: Standard error of the estimate
* `p_value`: Significance value of the estimate

## Implementation Details

### Required R Libraries
The script relies on the following R packages:

```R
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
```

### GLMM Model
The script applies a generalized linear mixed model (GLMM) with a binomial family:

```R
model <- glmer(refalt ~ offset(logit(ALT_dna_prob)) + toALT_edit_bias + 
                          (1 + toALT_edit_bias | Donor),
                          family = binomial, data = long_DF)
```

## Notes
* Ensure the input file is formatted correctly with tab-separated values
* The output results are saved as a TSV file at the specified `--output` path
* In this analysis, technical replicates were summed before processing

## License
This project is licensed under the MIT License.
