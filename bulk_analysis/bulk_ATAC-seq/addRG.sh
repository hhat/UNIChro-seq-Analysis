#!/bin/bash
#$ -S /bin/sh

export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

wd=${1}
sample_id=${2}


input_bam=${wd}/bowtie2/${sample_id}/${sample_id}_markdup_F1804f2q30_autosomal_bfilt.bam
out_bam=${wd}/bowtie2/${sample_id}/${sample_id}_markdup_F1804f2q30_autosomal_bfilt_addRG.bam

JAVA_TOOL_OPTIONS="-Xms12g -Xmx12g"

picard AddOrReplaceReadGroups \
    INPUT=${input_bam} \
    OUTPUT=${out_bam} \
    RGID=${sample_id} \
    RGLB=${sample_id} \
    RGPU=${sample_id} \
    RGPL=illumina \
    RGSM=${sample_id} \
    CREATE_INDEX=true

