#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
samplefile=info/2023-11-30.samples
mainodir=TMP/fastq/20231130
dd=/home/imgkono/data/img/miseq/20231130/Fastq
id=$( cat -n $samplefile | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
F1=$dd/${id}_R1_001.fastq.gz
F2=$dd/${id}_R2_001.fastq.gz

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $mainodir

# Extract UMI sequences (first 17 bases from Read1)
zcat $F1 |
awk '{if(NR%4==2){print }}' |
awk '{out = substr($1,1,17); print out}' > $mainodir/$id.umi

# Add UMI to read name (R1)
zcat $F1 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
paste $mainodir/$id.umi - |
awk 'BEGIN{FS="\t"}{
   umi=$1;
   readid=$2;
   sequence=$3;
   strand=$4;
   quality=$5;
   split(readid, D, " ");
   print D[1] "_" umi " " D[2] "\n" sequence "\n" strand "\n" quality
}' |
gzip -c - > $mainodir/$id.preqc.umi.R1.fastq.gz

# Add UMI to read name (R2)
zcat $F2 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
paste $mainodir/$id.umi - |
awk 'BEGIN{FS="\t"}{
   umi=$1;
   readid=$2;
   sequence=$3;
   strand=$4;
   quality=$5;
   split(readid, D, " ");
   print D[1] "_" umi " " D[2] "\n" sequence "\n" strand "\n" quality
}' |
gzip -c - > $mainodir/$id.preqc.umi.R2.fastq.gz

gzip -f $mainodir/$id.umi