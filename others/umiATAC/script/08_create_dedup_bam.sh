#!/bin/bash
#$ -S /bin/bash

sample_id=${1}
wd=${wd:-/home/ha7477/works/umi/umiATAC}
batch=${batch:-20251226}

if [ -z "$sample_id" ]; then
  echo "Usage: $0 <sample_id>"
  exit 1
fi

export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

bam=${wd}/bowtie/${batch}/${sample_id}/d1.sorted.q30.validumi.bam
odir=${wd}/umi_dedup/${sample_id}
mkdir -p "$odir"
obam=${odir}/${sample_id}_dedup.bam

samtools view "$bam" | awk '{
  split($1, a, "_")
  umi = a[2]
  chr = $3
  start = $4
  if(and($2, 16)) { strand = "-" } else { strand = "+" }
  key = chr":"start":"strand":"umi
  print key"\t"$1
}' | sort -k1,1 -S 8G | awk -F'\t' '
{
  if($1 != prev_key) {
    print $2
    prev_key = $1
  }
}' > "${odir}/${sample_id}_keep_reads.txt"

samtools view -h "$bam" | awk -v readfile="${odir}/${sample_id}_keep_reads.txt" '
BEGIN {
  while((getline line < readfile) > 0) {
    keep[line] = 1
  }
}
{
  if(/^@/) { print; next }
  if($1 in keep) { print }
}' | samtools view -bS - > "${obam}.tmp"

samtools sort "${obam}.tmp" -o "$obam"
rm "${obam}.tmp"
samtools index "$obam"

orig_reads=$(samtools view -c "$bam")
dedup_reads=$(samtools view -c "$obam")
echo "Sample: ${sample_id}"
echo "Original reads: ${orig_reads}"
echo "Dedup reads: ${dedup_reads}"
echo "Dedup ratio: $(echo "scale=2; ${dedup_reads}*100/${orig_reads}" | bc)%"

rm "${odir}/${sample_id}_keep_reads.txt"
