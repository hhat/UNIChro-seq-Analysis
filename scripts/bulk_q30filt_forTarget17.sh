#!/bin/bash
#$ -S /bin/bash

id=${1}

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH

IBAM=/home/ha7477/share/to_imgkono2/result/Jurkat/bowtie2/${id}/${id}_sorted.bam
OBAM=/home/ha7477/works/umi/Target17/bulk_reanalysis/${id}_sorted.q30.bam

samtools view -f 2 -q 30 $IBAM -o $OBAM
samtools index $OBAM
