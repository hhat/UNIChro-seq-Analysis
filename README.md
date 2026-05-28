# UNIChro-seq-Analysis

[![DOI](https://zenodo.org/badge/800069394.svg)](https://doi.org/10.5281/zenodo.15892171)

This repository contains the analysis scripts used in the UNIChro-seq
manuscript. Scripts are grouped by analysis module and numbered in the order in
which they were run. Raw sequencing data and large genome reference files are
not included; local paths, sample IDs, qsub array ranges, and reference files
should be edited before reuse.

## Folder Structure

* `Unichro-seq-pipeline/design_primer/`: primer design and in silico specificity checks.
* `Unichro-seq-pipeline/ATAC-seq/`: UNIChro-seq ATAC-seq preprocessing and allele-specific UMI counting.
* `Unichro-seq-pipeline/DNA-seq/`: DNA-seq preprocessing and allele-specific read counting.
* `crispr_analysis/CRISPResso2/`: CRISPR editing outcome analysis.
* `crispr_analysis/uni-directional_edit_demo/`: runnable demo for uni-directional editing analysis.
* `others/umiATAC/`: auxiliary UMI-ATAC comparison analysis.
* `others/siLEF1_RNA-seq/` and `others/siLEF1_ATAC-seq/`: auxiliary RNA-seq and ATAC-seq analysis scripts.

## Script Order

The analysis scripts have numbered file names (`01_`, `02_`, etc.) within each
analysis directory. These prefixes indicate the order in which the scripts were
run. `workflow.sh` files provide the corresponding command sequence for the main
UNIChro-seq primer-design, ATAC-seq, and DNA-seq workflows. Files beginning with
`00_` are helper wrappers. In `others/siLEF1_ATAC-seq/`, the `20_`-series files
are separate auxiliary scripts rather than a continuation of the main numbered
workflow.

## Software Requirements

The scripts were written for an SGE/qsub-based Linux HPC environment using
`bash`/`sh`, `Rscript`, conda-managed local environments, and standard Unix
tools. Recorded command-line tool versions are listed below.

| Tool | Version |
| --- | --- |
| bcl2fastq | v2.20.0.422 |
| FastQC | v0.11.9 |
| Cutadapt | v2.6 |
| Bowtie2 | v2.2.6 |
| Bowtie | v1.3.1 |
| Samtools | v1.9 |
| bedtools | v2.29.1 |
| Bcftools | v1.12 |
| Picard | v2.26.2 |
| MACS2 | v2.2.7.1 |
| GATK | v4.2.0.0 |
| STAR | v2.7.9a |
| RNA-SeQC2 | v2.4.2 |
| seqtk | v1.2-r94 |
| featureCounts | v2.0.3 |
| CRISPResso2 | v2.3.1 |
| tabix/bgzip | htslib v1.12 |

Recorded R package versions include:

| R package | Version |
| --- | --- |
| openPrimeR | v1.20.0 |
| DESeq2 | v1.34.0 |
| data.table | v1.16.4 |
| magrittr | v2.0.3 |
| dplyr | v1.1.4 |
| tidyr | v1.3.1 |
| tidyverse | v2.0.0 |
| ggplot2 | v3.5.1 |
| lme4 | v1.1-35.3 |
| lmerTest | v3.1-3 |
| ATACseqQC | v1.18.1 |
| Rsamtools | v2.14.0 |
| TxDb.Hsapiens.UCSC.hg38.knownGene | v3.16.0 |

R v4.2.1 and full package session information were recorded for the runnable
uni-directional editing demo in
`crispr_analysis/uni-directional_edit_demo/README.md`.

External reference files are required but not included in this repository,
including GRCh38 FASTA and Bowtie2 indexes, a STAR index, GENCODE gene
annotations, 1000 Genomes VCFs, a b37-to-hg38 chain file, ENCODE/hg38 blacklist
regions, and project-specific locus, probe, guide RNA, barcode, and sample
information files.

## License

This project is licensed under the terms of the LICENSE file included in this
repository.
