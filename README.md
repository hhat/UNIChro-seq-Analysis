# GLM Analysis Script

## Overview
This script processes a specified subset of columns from the input file and performs GLM analysis for each SNP using `ALT_dna_prob` as an offset. The output contains SNP-wise GLM results without `sample_base`, as per the updated requirement.

## Usage

### Running the Script
```bash
Rscript glm_analysis.R --input data/input_file.tsv --output results/glm_results.tsv
```

### Command-line Arguments
| Argument | Description |
|----------|-------------|
| `--input` | Path to the input TSV file |
| `--output` | Path to the output TSV file |

## Expected Input Format (`input_file.tsv`)
The script expects an input file with the following required columns:

| sample_id | SNP | ref | alt | REF_count | ALT_count | edit_direction |
|-----------|-----|-----|-----|------------|------------|----------------|
| R11_01_A | chr10_35126627_A_G | 248091 | 62899 | 60040789 | ALT_to_REF |
| R11_01_A | chr16_11088277_C_G | 235008 | 49465 | 12229536 | ALT_to_REF |
| R11_01_A | chr17_40096407_G_A | 79326 | 8691 | 03669693 | ALT_to_REF |
| R11_01_A | chr17_45895215_C_T | 42521 | 5688 | 917639938 | ALT_to_REF |
| R11_01_A | chr17_45895714_A_G | 132974 | 4416 | 542972277 | ALT_to_REF |
| R11_01_A | chr20_54173204_C_G | 193096 | 3111 | 09351338 | ALT_to_REF |

### Column Descriptions
* `sample_id`: Unique identifier for the sample
* `SNP`: SNP ID (chromosome and position)
* `ref`: Reference allele read count from ATAC-seq
* `alt`: Alternative allele read count from ATAC-seq
* `REF_count`: Reference allele read count from DNA-seq
* `ALT_count`: Alternative allele read count from DNA-seq
* `edit_direction`: Editing direction (`REF_to_ALT`, `ALT_to_REF`, or `NON_EDIT`)

## Output Format (`glm_results.tsv`)
The script generates an output file containing the results of the GLM analysis in the following format:

| SNP | effect | Estimate | Std_Error | p_value |
|-----|---------|-----------|------------|----------|
| chr10_35126627_A_G | Intercept | 0.1204 | 0.043 | 0.005 |
| chr10_35126627_A_G | toALT_edit_linear | -0.0852 | 0.024 | 0.001 |
| chr10_35126627_A_G | toALT_edit_nonlinear | 0.0501 | 0.020 | 0.012 |

### Output Column Descriptions
* `SNP`: The SNP ID being analyzed
* `effect`: The coefficient term in the GLM
  * `Intercept`: Baseline log odds ratio (how different ATAC-seq is from DNA-seq)
  * `toALT_edit_linear`: Linear effect of the edit direction
  * `toALT_edit_nonlinear`: Nonlinear effect of the edit direction
* `Estimate`: The estimated coefficient value
* `Std_Error`: Standard error of the estimate
* `p_value`: Significance value of the estimate

## Implementation Details

### Required R Libraries
The script relies on the following R packages:

```R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
```

### GLM Model
The script applies a generalized linear model (GLM) with a binomial family:

```R
model <- glm(refalt ~ offset(logit(ALT_count / (REF_count + ALT_count))) +
             toALT_edit_linear + toALT_edit_nonlinear,
             family = binomial, data = long_DF)
```

## Notes
* Ensure the input file is formatted correctly with tab-separated values (TSV)
* The script removes `sample_base` from the analysis as per the updated requirements
* The output results are saved as a TSV file at the specified `--output` path

## License
This project is licensed under the MIT License.
