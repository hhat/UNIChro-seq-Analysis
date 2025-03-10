# UNIChro-seq Bidirectional Analysis Demo

This demo-script identifies caQTLs and editing biases using a generalized linear model that combines UMI-counted ATAC-seq reads and DNA-seq ALT/(REF+ALT) ratios, both analyzed with allele-specific probes.

## Usage

### Running the Script
```bash
Rscript bidirectional_analysis.R --input data/input_file.tsv --output results/glm_results.tsv
```

### Command-line Arguments
| Argument | Description |
|----------|-------------|
| `--input` | Path to the input TSV file |
| `--output` | Path to the output TSV file |

## Expected Input Format (`input_file.tsv`)
The script expects an input file with the following required columns:

| Donor | SNP | ref | alt | REF_count | ALT_count | edit_direction |
|-----------|-----|-----|-----|------------|------------|----------------|
| R01 | rs111 | 24809 |10628 | 99600 | 40789 | ALT_to_REF |
| R01 | rs999 | 23500 | 8494 | 65122 | 29536 | ALT_to_REF |


### Column Descriptions
* `Donor`: Unique identifier for the Donor
* `SNP`: SNP ID (chromosome and position)
* `ref`: Reference allele read count from ATAC-seq
* `alt`: Alternative allele read count from ATAC-seq
* `REF_count`: Reference allele read count from DNA-seq
* `ALT_count`: Alternative allele read count from DNA-seq
* `edit_direction`: Editing direction (`REF_to_ALT`, `ALT_to_REF`, or `NON_EDIT`)

## Output Format (`results.tsv`)
The script generates an output file containing the results of the GLMM analysis in the following format:

| SNP | effect | Estimate | Std_Error | p_value |
|-----|---------|-----------|------------|----------|
| rs111 | Intercept | 0.1204 | 0.043 | 0.005 |
| rs111 | toALT_edit_linear | -0.0852 | 0.024 | 0.001 |

### Output Column Descriptions
* `SNP`: The SNP ID being analyzed
* `effect`: The coefficient term in the GLM
  * `Intercept`: Baseline log odds ratio (how different ATAC-seq is from DNA-seq)
  * `toALT_edit_linear`: Linear effect of the edit direction
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
model <- glmer(refalt ~ offset(logit(ALT_dna_prob)) + toALT_edit_linear + 
                          (1 + toALT_edit_linear | Donor),
                          family = binomial, data = long_DF)
```

## Notes
* Ensure the input file is formatted correctly with tab-separated values (TSV)
* The output results are saved as a TSV file at the specified `--output` path

## License
This project is licensed under the MIT License.
