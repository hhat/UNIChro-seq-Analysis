#!/bin/bash
#$ -S /bin/bash

fastq_dir=${1}
sample_file=${2}
output_dir=${3}

if [ -z "$fastq_dir" ] || [ -z "$sample_file" ] || [ -z "$output_dir" ]; then
  echo "Usage: $0 <bcl2fastq_fastq_dir> <sample_file> <output_dir>"
  exit 1
fi

mkdir -p "$output_dir"

while read sample_id; do
  echo "Merging lanes for $sample_id"
  cat "$fastq_dir"/${sample_id}_S*_L00*_R1_001.fastq.gz > "$output_dir"/${sample_id}_R1.fastq.gz
  cat "$fastq_dir"/${sample_id}_S*_L00*_R2_001.fastq.gz > "$output_dir"/${sample_id}_R2.fastq.gz
done < "$sample_file"
