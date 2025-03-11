#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
sample_barcode_file=info/2024-06-17.samples
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
dd=/home/imgkono/data/img/miseq/20240617/Fastq
ODIR=cutadapt/20240617/$id

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $ODIR

# Process reads with Cutadapt
IF1=$dd/${id}_R1_001.fastq.gz
IF2=$dd/${id}_R2_001.fastq.gz
OF1=$ODIR/$id.postqc.R1.fastq.gz
OF2=$ODIR/$id.postqc.R2.fastq.gz

# Trim adapter sequences
cutadapt \
   -a CTGTCTCTTATACACATCT \
   -A CAGAGATCGGAAGAGCG \
   --minimum-length 20 \
   -o  $OF1 -p  $OF2 \
   $IF1 $IF2

# Count reads before and after trimming
N_Input_R1=$( zcat $IF1 | wc -l | awk '{ print $1 /4 }' )
N_Input_R2=$( zcat $IF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )

# Save read count summary
echo $id $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 >> $ODIR/round1.count