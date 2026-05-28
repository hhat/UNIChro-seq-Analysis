#!/bin/sh

jobid=${SGE_TASK_ID}

wd=${wd:-/home/ha7477/works/umi/umiATAC}
batch=${batch:-20251226}
samplefile=${samplefile:-${wd}/info/2025-12-26.samples}

sample_id=$( cat -n "$samplefile" | awk -v jobid="$jobid" '{if($1==jobid){print $2}}' )
IBAM=${wd}/bowtie/${batch}/${sample_id}/d1.sorted.q30.bam
OBAM=${wd}/bowtie/${batch}/${sample_id}/d1.sorted.q30.validumi.bam

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

samtools view -h "$IBAM" | grep -v "_NA" | samtools view -bS - > "$OBAM"
samtools index "$OBAM"

N_before=$( samtools view -c "$IBAM" )
N_after=$( samtools view -c "$OBAM" )
N_filtered=$( expr "$N_before" - "$N_after" )

echo "$sample_id $N_before $N_after $N_filtered" >> "${wd}/bowtie/${batch}/na_umi_filter_stats.txt"
