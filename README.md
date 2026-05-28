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
* `others/siLEF1_RNA-seq/` and `others/siLEF1_ATAC-seq/`: auxiliary siLEF1 RNA-seq and ATAC-seq analyses.

## Script Order

The analysis scripts have numbered file names (`01_`, `02_`, etc.) within each
analysis directory. These prefixes indicate the order in which the scripts were
run. `workflow.sh` files provide the corresponding command sequence for the main
UNIChro-seq primer-design, ATAC-seq, and DNA-seq workflows. Files beginning with
`00_` are helper wrappers, and `20_`-series siLEF1 ATAC-seq scripts are
Target17/figure-support analyses.

## Software Requirements

The scripts were written for an SGE/qsub-based Linux HPC environment using
`bash`/`sh`, `Rscript`, and standard Unix tools. The main command-line tools are
`bcl2fastq` 2.20.0.422, FastQC, Cutadapt v2.6, Bowtie2 v2.2.6, Bowtie v1.3.1,
Samtools v1.9, bedtools, Bcftools v1.12, Picard v2.26.2, MACS2 v2.2.7.1, GATK
4.2.0.0, STAR v2.7.9a, RNA-SeQC2 v2.4.2, seqtk, featureCounts, CRISPResso2,
tabix/bgzip, and conda-managed local environments.

R packages used by the scripts include openPrimeR v1.20.0, DESeq2 v1.34.0,
data.table, magrittr, dplyr, tidyr, tidyverse, ggplot2, lme4, lmerTest,
ATACseqQC, Rsamtools, and TxDb.Hsapiens.UCSC.hg38.knownGene. R 4.2.1 and the
package versions listed in `crispr_analysis/uni-directional_edit_demo/README.md`
were recorded for the runnable uni-directional editing demo.

External reference files are required but not included in this repository,
including GRCh38 FASTA and Bowtie2 indexes, a STAR index, GENCODE gene
annotations, 1000 Genomes VCFs, a b37-to-hg38 chain file, ENCODE/hg38 blacklist
regions, and project-specific locus, probe, guide RNA, barcode, and sample
information files.

## License

This project is licensed under the terms of the LICENSE file included in this
repository.
