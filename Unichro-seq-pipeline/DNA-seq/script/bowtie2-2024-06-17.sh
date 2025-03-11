#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
threads=4  # Number of alignment threads to launch
sample_barcode_file=info/2024-06-17.samples
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
F1=cutadapt/20240617/$id/$id.postqc.R1.fastq.gz
F2=cutadapt/20240617/$id/$id.postqc.R2.fastq.gz
ODIR=bowtie/20240617/$id
OFILE=bowtie/20240617/$id/d1
Index=/home/imgishi/reference/Bowtie2_ref/GRCh38/GRCh38
IBAM=bowtie/20240617/$id/d1.sorted.bam
OBAM=bowtie/20240617/$id/d1.sorted.q30.bam

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

# Process and filter alignment results
samtools view -bS $OFILE.sam > $OFILE.bam

# Filter for properly paired reads with MAPQ >= 10
samtools view -f 2 -q 10 $OFILE.bam -o $OFILE.qc.bam

samtools sort $OFILE.qc.bam -o $OFILE.sorted.bam
samtools index $OFILE.sorted.bam

# Apply MAPQ 30 filtering
samtools view -f 2 -q 30 $IBAM -o $OBAM
samtools index $OBAM

# Remove intermediate files
rm -f $OFILE.sam
rm -f $OFILE.qc.bam