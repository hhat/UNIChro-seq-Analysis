#!/bin/bash
#$ -S /bin/bash

wd=$1
sample_id=$2

fastq1=${wd}/merged_fastq/${sample_id}_R1.fastq.gz
fastq2=${wd}/merged_fastq/${sample_id}_R2.fastq.gz

export PATH=/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

mkdir -p ${wd}/trim_fastq/
trim_fastq1=${wd}/trim_fastq/${sample_id}_val_1.fq.gz
trim_fastq2=${wd}/trim_fastq/${sample_id}_val_2.fq.gz

cpu=2
cutadapt -j ${cpu} -e 0.1 -q 0 --minimum-length 20 --pair-filter=any -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${trim_fastq1} -p ${trim_fastq2} ${fastq1} ${fastq2}

