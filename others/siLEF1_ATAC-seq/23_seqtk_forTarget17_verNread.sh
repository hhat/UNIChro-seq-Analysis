#!/bin/bash
#$ -S /bin/bash

wd=${1}
sample_id=${2}
num_reads=${3}

seed=$(echo "1234 + ($num_reads * 1000)" | bc)
seed=${seed%.*} 

export PATH=/home/ha7477/tools/miniconda3/envs/de/bin::${PATH}

R1=${wd}/trim_fastq/${sample_id}_val_1.fq.gz
R2=${wd}/trim_fastq/${sample_id}_val_2.fq.gz
od=${wd}/downsampling_fastq
mkdir -p ${od}/${num_reads}

/home/ha7477/tools/miniconda3/envs/de/bin/seqtk sample -s${seed} ${R1} ${num_reads} | gzip > ${od}/${num_reads}/${sample_id}_${num_reads}_R1.fastq.gz
/home/ha7477/tools/miniconda3/envs/de/bin/seqtk sample -s${seed} ${R2} ${num_reads} | gzip > ${od}/${num_reads}/${sample_id}_${num_reads}_R2.fastq.gz


