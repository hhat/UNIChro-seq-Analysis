#!/bin/bash
#$ -S /bin/sh


export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

wd=${1}
sample_id=${2}
core=8

bam=${wd}/bowtie2/${sample_id}/${sample_id}_sorted.bam
outbam=/home/ha7477/share/to_imgkono2/result/Jurkat/Target17/${sample_id}_f2q30_sorted.bam

##samtools 
samtools view -f 2 -q 30 ${bam} -Obam | samtools sort - -o ${outbam}
samtools index ${outbam}
