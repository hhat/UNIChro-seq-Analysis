# UNIChro-seq Nested Multiplex Primer Design

This repository provides a pipeline for designing nested multiplex primers for UNIChro-seq.

## Overview

This pipeline designs optimal nested multiplex primers for specified SNPs and evaluates their mapping characteristics. Primer design is performed using the [openPrimeR](https://www.bioconductor.org/packages/release/bioc/vignettes/openPrimeR/inst/doc/openPrimeR_vignette.html) package.

## Required Libraries

### R packages:
- openPrimeR (R version 4.3.3) also requiring MELTING, OligoArrayAux

### Command-line tools:
- samtools
- bowtie
- bedtools
- perl


## Usage

### Input File Format

Example of `demo.input`:
```
chr1_12040655_C_T chr1_12040288_12040958 1
chr1_183186264_A_G chr1_183185574_183186564 1
chr1_202174197_G_C chr1_202174106_202175063 2
chr1_209611861_G_A chr1_209611564_209612178 2
```

Column descriptions:
- Column 1: SNP information (chr_pos_ref_alt)
- Column 2: Open chromatin region (chr_from_to)
- Column 3: Direction (1=forward strand, 2=reverse strand)
- Column 4: Distance from SNP (Dist_snp parameter)


```bash
./workflow/run_full_pipeline.sh example/demo.input results/my_primers "Tm=62;Window=50;CVG=0.8"
```

Alternatively, you can run each step individually:

### STEP 1: Prepare Input Data

Prepare your input file:
```bash
mkdir -p results/my_primers
cp example/demo.input results/my_primers/
```

### STEP 2: Create FASTA Files

```bash
script/prepare_fastq.sh results/my_primers/demo.input results/my_primers
```

This step extracts 200bp genomic sequences containing each SNP position and generates the `target.200bp.fa` file.

### STEP 3: Design Primers

```bash
script/design_primer.sh "Tm=62;Window=50;CVG=0.8;ODIR=results/my_primers"
```

Parameter descriptions:
- `Tm`: Target melting temperature for primers (±3°C)
- `Window`: Size of the region for primer design
- `CVG`: Stringency of design (0-1). 1 means design primers for all SNPs by relaxing constraints, 0 means maintain strict conditions
- `ODIR`: Output directory

### STEP 4: Evaluate Primer Mapping

```bash
script/evaluate_primer.sh results/my_primers
```

This step evaluates the mapping characteristics of the designed primers to the genome.

### STEP 5: Summarize Results

```bash
Rscript script/make_summary.R results/my_primers/demo.input results/my_primers
```

This step integrates all information and generates the final result files `target.primers.mapinfo.txt` and `target.primers.mapinfo_oneline.txt`.

## Output Files

- `target.200bp.fa`: 200bp genomic sequences containing each SNP position
- `target.primers`: Information about designed primers
- `g38_map.allseq.summary.txt`: Mapping information for the full primer sequences
- `g38_map.15.summary.txt`: Mapping information for the 3' end 15bp of primers
- `target.primers.mapinfo.txt`: Integrated result file with all information
- `target.primers.mapinfo_oneline.txt`: Result file with one line per SNP

## Parameter Optimization

The most critical parameter is `CVG`. This value determines the balance between primer design stringency and success rate:

- `CVG=1.0`: Relaxes conditions to design primers for all SNPs
- `CVG=0.0`: Maintains strict conditions and designs primers only for SNPs that meet the criteria

We recommend adjusting the `CVG` value based on the characteristics of your SNPs and project requirements.

## References
Kreer C, Döring M, Lehnen N, Ercanoglu M, Gieselmann L, Luca D, Jain K, Schommers P, Pfeifer N, Klein F (2020). “openPrimeR for multiplex amplification of highly diverse templates.” Journal of Immunological Methods, 112752. https://bioconductor.org/packages/release/bioc/html/openPrimeR.html.

Döring M, Kreer C, Lehnen N, Klein F, Pfeifer N (2019). “Modeling the Amplification of Immunoglobulins through Machine Learning on Sequence-Specific Features.” Scientific Reports, 9(1), 1–11.
