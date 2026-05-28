#!/bin/bash
#$ -S /bin/bash

sample_id=${1}
wd=${wd:-/home/ha7477/works/umi/umiATAC}
core=${core:-8}
blacklist=${blacklist:-/home/ha7477/reference/atac/ATACblacklist/hg38.blacklist.bed.gz}

if [ -z "$sample_id" ]; then
  echo "Usage: $0 <sample_id>"
  exit 1
fi

export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

inputbam=${wd}/umi_dedup/${sample_id}/${sample_id}_dedup.bam
odir=${wd}/umi_dedup_filtered/${sample_id}
mkdir -p "$odir"

samtools flagstat "$inputbam" > "${odir}/${sample_id}_dedup.bam_flagstat"
samtools view -F 780 -f 2 -q 30 "$inputbam" -Obam | \
  samtools sort -@ "$core" -o "${odir}/${sample_id}_F780f2q30.bam"

samtools flagstat "${odir}/${sample_id}_F780f2q30.bam" > "${odir}/${sample_id}_F780f2q30.bam_flagstat"
samtools index "${odir}/${sample_id}_F780f2q30.bam"

outputbam=${odir}/${sample_id}_F780f2q30_autosomal_bfilt.bam

samtools view -o - "${odir}/${sample_id}_F780f2q30.bam" $(seq 1 22 | sed 's/^/chr/') -b | \
  bedtools intersect -nonamecheck -v -abam "stdin" -b "$blacklist" | \
  samtools sort -@ "$core" > "$outputbam"

samtools index "$outputbam"
samtools flagstat "$outputbam" > "${outputbam}_flagstat"

rm "${odir}/${sample_id}_F780f2q30.bam"
rm "${odir}/${sample_id}_F780f2q30.bam.bai"

echo "Done: ${sample_id}"
echo "Final output: ${outputbam}"
