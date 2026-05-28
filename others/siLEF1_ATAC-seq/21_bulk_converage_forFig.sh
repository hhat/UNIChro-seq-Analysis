#!/bin/bash
#$ -S /bin/sh

wd=${1}
sample_id=${2}

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
  
bed="/home/ha7477/works/umi/Target17/genome_coverage/grch38_1e7_bin.txt"

#paird R2 read=130
BAM=/home/ha7477/share/to_imgkono2/result/Jurkat/Target17/${sample_id}_f2q30_sorted.bam
samtools view -@ 8 -b -f 130 $BAM chr{1..22} |
   bedtools coverage -a ${bed} -b - > /home/ha7477/share/to_imgkono2/result/Jurkat/Target17/genome_coverage/${sample_id}_coverage.txt
