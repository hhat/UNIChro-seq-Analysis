#!/bin/bash
#$ -S /bin/bash

wd=${1}
sample_id=${2}
cpu=16

export PATH=/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

R1=${wd}/trim_fastq/${sample_id}_val_1.fq.gz
R2=${wd}/trim_fastq/${sample_id}_val_2.fq.gz

Bowtie2Index=/home/ha7477/reference/atac/bowtie2/references_Homo_sapiens_assembly38_noALT_noHLA_noDecoy


tmpoutdir=${wd}/bowtie2/${sample_id}

mkdir -p $tmpoutdir
cd $tmpoutdir

bowtie2 -k 4 -X 2000 --sensitive --local --no-unal -p $cpu -x $Bowtie2Index -1 ${R1} -2 ${R2} 2> ${sample_id}.bowtie2.error | 
samtools sort -@ $cpu -O bam -o ${tmpoutdir}/${sample_id}_sorted.bam

samtools index -@ $cpu ${tmpoutdir}/${sample_id}_sorted.bam

