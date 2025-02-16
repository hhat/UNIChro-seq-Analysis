#!/bin/bash
#$ -S /bin/bash

wd=${1}
sample_id=${2}
num_reads=${3}
base_num_reads=${4}

seed=$(echo "1234 + ($num_reads * 1000)" | bc)
seed=${seed%.*}  # 小数点以下を削除

export PATH=/home/ha7477/tools/miniconda3/envs/de/bin::${PATH}

od=${wd}/downsampling_fastq
mkdir -p ${od}/${num_reads}
R1=${od}/${base_num_reads}/${sample_id}_${base_num_reads}_R1.fastq.gz
R2=${od}/${base_num_reads}/${sample_id}_${base_num_reads}_R2.fastq.gz


# num_readsで指定したリード数にダウンサンプリング
/home/ha7477/tools/miniconda3/envs/de/bin/seqtk sample -s${seed} ${R1} ${num_reads} | gzip > ${od}/${num_reads}/${sample_id}_${num_reads}_R1.fastq.gz
/home/ha7477/tools/miniconda3/envs/de/bin/seqtk sample -s${seed} ${R2} ${num_reads} | gzip > ${od}/${num_reads}/${sample_id}_${num_reads}_R2.fastq.gz

