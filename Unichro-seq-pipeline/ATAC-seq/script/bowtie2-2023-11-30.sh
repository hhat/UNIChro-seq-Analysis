#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
threads=4  # Number of alignment threads to launch
sample_barcode_file=info/2023-11-30.sample_barcode_pair
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
F1=cutadapt/20231130/$id/$id.$barcode.postqc.umi.R1.fastq.gz
F2=cutadapt/20231130/$id/$id.$barcode.postqc.umi.R2.fastq.gz
ODIR=bowtie/20231130/$id.$barcode
OFILE=bowtie/20231130/$id.$barcode/d1
Index=/home/imgishi/reference/Bowtie2_ref/GRCh38/GRCh38

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $ODIR

# alignment with Bowtie2
bowtie2 -x $Index \
   -1 $F1 \
   -2 $F2 \
   --threads $threads \
   --no-discordant \
   --maxins 1000 \
   -S $OFILE.sam

# filter mapping output
samtools view -bS -f 2 -q 30 $OFILE.sam > $OFILE.q30.bam

samtools sort $OFILE.q30.bam -o $OFILE.sorted.q30.bam
samtools index $OFILE.sorted.q30.bam

# Remove intermediate files
rm -f $OFILE.sam
rm -f $OFILE.q30.bam