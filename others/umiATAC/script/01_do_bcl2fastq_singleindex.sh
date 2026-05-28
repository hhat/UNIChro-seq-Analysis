#!/bin/bash
#$ -S /bin/bash

runfolder_dir=${1}
sample_sheet=${2}
output_dir=${3}

if [ -z "$runfolder_dir" ] || [ -z "$sample_sheet" ] || [ -z "$output_dir" ]; then
  echo "Usage: $0 <runfolder_dir> <sample_sheet.csv> <output_dir>"
  exit 1
fi

/usr/local/package/bcl2fastq/2.20.0.422/bin/bcl2fastq \
  --runfolder-dir "$runfolder_dir" \
  --sample-sheet "$sample_sheet" \
  --output-dir "$output_dir" \
  --use-bases-mask Y100n,I8,Y100n \
  --processing-threads 28 \
  --loading-threads 4 \
  --writing-threads 4
