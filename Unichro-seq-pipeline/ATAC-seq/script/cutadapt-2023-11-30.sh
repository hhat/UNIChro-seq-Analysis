#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
sample_barcode_file=info/2023-11-30.sample_barcode_pair
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
IDIR=demultiplex/20231130/$id
ODIR=cutadapt/20231130/$id

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $ODIR

# Round 1: Trim AGATGTGTATAAGAGACAG from 5' ONLY ONCE
# Outputs reads without AGATGTGTATAAGAGACAG to diagnose the experiments
IF1=$IDIR/$id.$barcode.preqc.umi.R1.fastq.gz
IF2=$IDIR/$id.$barcode.preqc.umi.R2.fastq.gz
OF1=$ODIR/$id.$barcode.postqc_tmpA.umi.R1.fastq.gz
OF2=$ODIR/$id.$barcode.postqc_tmpA.umi.R2.fastq.gz
OF1_notrim=$ODIR/$id.$barcode.postqc_tmpA_notrim.umi.R1.fastq.gz
OF2_notrim=$ODIR/$id.$barcode.postqc_tmpA_notrim.umi.R2.fastq.gz

cutadapt \
   -g "AGATGTGTATAAGAGACAG;min_overlap=19" \
   --times 1 \
   --untrimmed-output  $OF1_notrim \
   --untrimmed-paired-output  $OF2_notrim \
   -o  $OF1 -p  $OF2 \
   $IF1 $IF2

N_Input_R1=$( zcat $IF1 | wc -l | awk '{ print $1 /4 }' )
N_Input_R2=$( zcat $IF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R1=$( zcat $OF1_notrim | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R2=$( zcat $OF2_notrim | wc -l | awk '{ print $1 /4 }' )

echo $id $barcode $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 $N_Output_notrim_R1 $N_Output_notrim_R2 >> $ODIR/round1.count

# Round 2: Trim AGATGTGTATAAGAGACAG from 5' MULTIPLE TIMES
# Exclude reads with multiple ME sequences
# In this round, the reads WITHOUT trimming are the ones we should use in downstream analyses
IF1=$ODIR/$id.$barcode.postqc_tmpA.umi.R1.fastq.gz
IF2=$ODIR/$id.$barcode.postqc_tmpA.umi.R2.fastq.gz
OF1=$ODIR/$id.$barcode.postqc_tmpB.umi.R1.fastq.gz
OF2=$ODIR/$id.$barcode.postqc_tmpB.umi.R2.fastq.gz
OF1_notrim=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R1.fastq.gz
OF2_notrim=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R2.fastq.gz

cutadapt \
   -g "AGATGTGTATAAGAGACAG;min_overlap=19" \
   --times 5 \
   --untrimmed-output  $OF1_notrim \
   --untrimmed-paired-output  $OF2_notrim \
   -o  $OF1 -p  $OF2 \
   $IF1 $IF2

N_Input_R1=$( zcat $IF1 | wc -l | awk '{ print $1 /4 }' )
N_Input_R2=$( zcat $IF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R1=$( zcat $OF1_notrim | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R2=$( zcat $OF2_notrim | wc -l | awk '{ print $1 /4 }' )

echo $id $barcode $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 $N_Output_notrim_R1 $N_Output_notrim_R2 >> $ODIR/round2.count

# Round 3: Trim CTGTCTCTTATACACATCT from 3' MULTIPLE TIMES
# Exclude reads with multiple reverse-complement ME sequences
IF1=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R1.fastq.gz
IF2=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R2.fastq.gz
OF1=$ODIR/$id.$barcode.postqc.umi.R1.fastq.gz
OF2=$ODIR/$id.$barcode.postqc.umi.R2.fastq.gz

cutadapt \
   -a CTGTCTCTTATACACATCT \
   -A CTGTCTCTTATACACATCT \
   --minimum-length 20 \
   --times 5 \
   -o  $OF1 -p  $OF2 \
   $IF1 $IF2

N_Input_R1=$( zcat $IF1 | wc -l | awk '{ print $1 /4 }' )
N_Input_R2=$( zcat $IF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )

echo $id $barcode $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 > $ODIR/round3.count

# Adapter sequences:
 #AGATGTGTATAAGAGACAG (R1,5'): -g (ME in a positive strand; in a standard library, this is not sequence because this works as a sequence primer, but in our pipeline, we do not use illumina primer cocktail)
# CTGTCTCTTATACACATCT (R1,3'): -a (ME in a reverse strand)
# CTGTCTCTTATACACATCT (R2,3'): -A (ME in a reverse strand)