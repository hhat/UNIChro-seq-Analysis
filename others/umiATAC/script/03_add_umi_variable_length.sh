#!/bin/sh

jobid=${SGE_TASK_ID}

wd=${wd:-/home/ha7477/works/umi/umiATAC}
batch=${batch:-20251226}
samplefile=${samplefile:-${wd}/info/2025-12-26.samples}
input_dir=${input_dir:-${wd}/fastq_merged}
output_dir=${output_dir:-${wd}/TMP/fastq/${batch}}
me_sequence=${me_sequence:-AGATGTGTATAAGAGACAG}

sample_id=$( cat -n "$samplefile" | awk -v jobid="$jobid" '{if($1==jobid){print $2}}' )
F1=${input_dir}/${sample_id}_R1.fastq.gz
F2=${input_dir}/${sample_id}_R2.fastq.gz

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p "$output_dir"

# Extract the variable-length UMI before the first ME sequence in Read 1.
# Valid UMI lengths are 6, 13, 20, and 25 bp. Other reads are marked as NA.
zcat "$F1" |
awk -v me="$me_sequence" '{
  if(NR%4==2){
    pos = index($0, me)
    if(pos > 0){
      umi_len = pos - 1
      if(umi_len == 6 || umi_len == 13 || umi_len == 20 || umi_len == 25){
        print substr($0, 1, pos - 1)
      } else {
        print "NA"
      }
    } else {
      print "NA"
    }
  }
}' > "$output_dir"/${sample_id}.umi

zcat "$F1" |
awk '{
  printf("%s",$0);
  if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
paste "$output_dir"/${sample_id}.umi - |
awk 'BEGIN{FS="\t"}{
  umi=$1;
  readid=$2;
  sequence=$3;
  strand=$4;
  quality=$5;
  split(readid, D, " ");
  print D[1] "_" umi " " D[2] "\n" sequence "\n" strand "\n" quality
}' |
gzip -c - > "$output_dir"/${sample_id}.preqc.umi.R1.fastq.gz

zcat "$F2" |
awk '{
  printf("%s",$0);
  if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
paste "$output_dir"/${sample_id}.umi - |
awk 'BEGIN{FS="\t"}{
  umi=$1;
  readid=$2;
  sequence=$3;
  strand=$4;
  quality=$5;
  split(readid, D, " ");
  print D[1] "_" umi " " D[2] "\n" sequence "\n" strand "\n" quality
}' |
gzip -c - > "$output_dir"/${sample_id}.preqc.umi.R2.fastq.gz

gzip -f "$output_dir"/${sample_id}.umi
