#!/bin/bash
#$ -S /bin/bash

sample_id=$1
R1=$2
R2=$3
resultdir=$4


export PATH=/home/ha7477/tools/miniconda3/envs/de/bin:/home/ha7477/tools/miniconda3/bin:${PATH}


##FastQC needs premaking of outputdir
mkdir -p ${resultdir}/R1
mkdir -p ${resultdir}/R2


fastqc --nogroup --outdir ${resultdir}/R1 ${R1}
fastqc --nogroup --outdir ${resultdir}/R2 ${R2}
