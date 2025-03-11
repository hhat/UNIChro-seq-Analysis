# UNIChro-seq Bidirectional Analysis Demo

This demo-script identifies caQTLs and editing biases using a generalized linear mixed model that combines UMI-counted ATAC-seq read counts and DNA-seq ALT/(REF+ALT) ratios, both analyzed with allele-specific probes.

## Usage

### Running the Script
```bash
Rscript bidirectional_analysis.R --input input_file.txt --output results.txt
```

### Command-line Arguments
| Argument | Description |
|----------|-------------|
| `--input` | Path to the input text file |
| `--output` | Path to the output text file |

## Expected Input Format (`input_file.txt`)
The script expects an tab-separated input file with the following required columns:

| Donor | edit_direction | SNP | ref | alt | REF_count | ALT_count |
|--------|---------------|-----|-----|-----|-----------|-----------| 
| S01 | ALT_to_REF | rs2248137 | 1930 | 963 | 111093 | 51338 |
| S01 | NON_EDIT | rs2248137 | 1303 | 1423 | 76226 | 57463 |
| S01 | REF_to_ALT | rs2248137 | 1512 | 2016 | 68108 | 74336 |
| S02 | ALT_to_REF | rs2248137 | 917 | 527 | 127240 | 55054 |
| S02 | NON_EDIT | rs2248137 | 1016 | 908 | 99598 | 83489 |
| S02 | REF_to_ALT | rs2248137 | 900 | 1188 | 65335 | 71255 |
| S03 | ALT_to_REF | rs2248137 | 1002 | 707 | 90922 | 46997 |

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
| rs2248137 | caQTL | 0.201 | 0.0528 | 0.000144 |
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
* The output results are saved as a tab-separated text file at the specified `--output` path
* In this analysis, technical replicates were summed before processing

## License
This project is licensed under the MIT License.
