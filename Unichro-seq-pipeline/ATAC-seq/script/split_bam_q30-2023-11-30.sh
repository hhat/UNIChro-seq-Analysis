#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
sample_barcode_file=info/2023-11-30.sample_barcode_pair
locus_info=info/2023-11-30.locus_info2
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
ODIR=bowtie/20231130/$id.$barcode
BAM=bowtie/20231130/$id.$barcode/d1.sorted.q30.bam

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1


cat $locus_info | awk '{print $1}' | while read snpid; do
   # Extract chromosome and position information
   CHR=$( echo $snpid | cut -d "_" -f 1 )
   POS=$( echo $snpid | cut -d "_" -f 2 )
   
   # Define region of interest (±100bp from SNP position)
   FROM=$( expr $POS - 100 )
   TO=$( expr $POS + 100 )
   if [ $FROM -lt 1 ] ; then FROM=1 ; fi
   Region=$CHR:${FROM}-${TO}
   
   mkdir -p $ODIR/$snpid
   
   # Extract reads from the region
   samtools view -b $BAM $Region > $ODIR/$snpid/$snpid.100bp_q30.bam

   samtools index $ODIR/$snpid/$snpid.100bp_q30.bam
done