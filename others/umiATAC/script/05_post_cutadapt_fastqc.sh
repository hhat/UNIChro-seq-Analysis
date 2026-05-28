#!/bin/sh

jobid=${SGE_TASK_ID}

wd=${wd:-/home/ha7477/works/umi/umiATAC}
batch=${batch:-20251226}
samplefile=${samplefile:-${wd}/info/2025-12-26.samples}
mainodir=${mainodir:-${wd}/fastqc/${batch}}

sample_id=$( cat -n "$samplefile" | awk -v jobid="$jobid" '{if($1==jobid){print $2}}' )
IDIR=${wd}/cutadapt/${batch}/${sample_id}
F1=${IDIR}/${sample_id}.postqc.umi.R1.fastq.gz
F2=${IDIR}/${sample_id}.postqc.umi.R2.fastq.gz

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

odir=${mainodir}/${sample_id}.R1.postqc
mkdir -p "$odir"
fastqc -o "$odir" -f fastq "$F1"

odir=${mainodir}/${sample_id}.R2.postqc
mkdir -p "$odir"
fastqc -o "$odir" -f fastq "$F2"
