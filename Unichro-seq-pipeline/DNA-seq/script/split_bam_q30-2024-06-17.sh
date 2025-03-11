#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
samplefile=info/2024-06-17.samples
id=$( cat -n $samplefile | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
locus_info=info/R10.locus_info2
ODIR=bowtie/20240617/$id
BAM=bowtie/20240617/$id/d1.sorted.q30.bam

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1


cat $locus_info | awk '{print $1}' | while read snpid; do
   # Extract chromosome and position information
   CHR=$( echo $snpid | cut -d "_" -f 1 )
   POS=$( echo $snpid | cut -d "_" -f 2 )
   
   # Define region of interest (±300bp from SNP position)
   FROM=$( expr $POS - 300 )
   TO=$( expr $POS + 300 )
   if [ $FROM -lt 1 ]; then FROM=1; fi
   Region=$CHR:${FROM}-${TO}
   
   mkdir -p $ODIR/$snpid
   
   # Extract reads from the region
   samtools view -b $BAM $Region > $ODIR/$snpid/$snpid.100bp_q30.bam
   
   samtools index $ODIR/$snpid/$snpid.100bp_q30.bam
done