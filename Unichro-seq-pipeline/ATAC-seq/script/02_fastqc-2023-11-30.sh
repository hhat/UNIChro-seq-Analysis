#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
samplefile=info/2023-11-30.samples
mainodir=fastqc/20231130
dd=/home/imgkono/data/img/miseq/20231130/Fastq
id=$( cat -n $samplefile | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
F1=$dd/${id}_R1_001.fastq.gz
F2=$dd/${id}_R2_001.fastq.gz

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

# Process R1 reads
odir=$mainodir/$id.R1
mkdir -p $odir
fastqc -o $odir -f fastq $F1

# Process R2 reads
odir=$mainodir/$id.R2
mkdir -p $odir
fastqc -o $odir -f fastq $F2