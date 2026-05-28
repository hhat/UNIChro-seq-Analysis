#!/bin/sh

jobid=${SGE_TASK_ID}

wd=${wd:-/home/ha7477/works/umi/umiATAC}
batch=${batch:-20251226}
samplefile=${samplefile:-${wd}/info/2025-12-26.samples}
input_dir=${input_dir:-${wd}/TMP/fastq/${batch}}
output_root=${output_root:-${wd}/cutadapt/${batch}}

sample_id=$( cat -n "$samplefile" | awk -v jobid="$jobid" '{if($1==jobid){print $2}}' )
ODIR=${output_root}/${sample_id}

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p "$ODIR"

IF1=${input_dir}/${sample_id}.preqc.umi.R1.fastq.gz
IF2=${input_dir}/${sample_id}.preqc.umi.R2.fastq.gz
OF1=${ODIR}/${sample_id}.postqc_tmpA.umi.R1.fastq.gz
OF2=${ODIR}/${sample_id}.postqc_tmpA.umi.R2.fastq.gz
OF1_notrim=${ODIR}/${sample_id}.postqc_tmpA_notrim.umi.R1.fastq.gz
OF2_notrim=${ODIR}/${sample_id}.postqc_tmpA_notrim.umi.R2.fastq.gz

cutadapt \
  -g "AGATGTGTATAAGAGACAG;min_overlap=19" \
  --times 1 \
  --untrimmed-output "$OF1_notrim" \
  --untrimmed-paired-output "$OF2_notrim" \
  -o "$OF1" -p "$OF2" \
  "$IF1" "$IF2"

echo "$sample_id" \
  "$( zcat "$IF1" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$IF2" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF1" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF2" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF1_notrim" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF2_notrim" | wc -l | awk '{ print $1 /4 }' )" >> "$ODIR"/round1.count

IF1=${ODIR}/${sample_id}.postqc_tmpA.umi.R1.fastq.gz
IF2=${ODIR}/${sample_id}.postqc_tmpA.umi.R2.fastq.gz
OF1=${ODIR}/${sample_id}.postqc_tmpB.umi.R1.fastq.gz
OF2=${ODIR}/${sample_id}.postqc_tmpB.umi.R2.fastq.gz
OF1_notrim=${ODIR}/${sample_id}.postqc_tmpB_notrim.umi.R1.fastq.gz
OF2_notrim=${ODIR}/${sample_id}.postqc_tmpB_notrim.umi.R2.fastq.gz

cutadapt \
  -g "AGATGTGTATAAGAGACAG;min_overlap=19" \
  --times 5 \
  --untrimmed-output "$OF1_notrim" \
  --untrimmed-paired-output "$OF2_notrim" \
  -o "$OF1" -p "$OF2" \
  "$IF1" "$IF2"

echo "$sample_id" \
  "$( zcat "$IF1" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$IF2" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF1" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF2" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF1_notrim" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF2_notrim" | wc -l | awk '{ print $1 /4 }' )" >> "$ODIR"/round2.count

IF1=${ODIR}/${sample_id}.postqc_tmpB_notrim.umi.R1.fastq.gz
IF2=${ODIR}/${sample_id}.postqc_tmpB_notrim.umi.R2.fastq.gz
OF1=${ODIR}/${sample_id}.postqc.umi.R1.fastq.gz
OF2=${ODIR}/${sample_id}.postqc.umi.R2.fastq.gz

cutadapt \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  --minimum-length 20 \
  --times 5 \
  -o "$OF1" -p "$OF2" \
  "$IF1" "$IF2"

echo "$sample_id" \
  "$( zcat "$IF1" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$IF2" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF1" | wc -l | awk '{ print $1 /4 }' )" \
  "$( zcat "$OF2" | wc -l | awk '{ print $1 /4 }' )" >> "$ODIR"/round3.count
