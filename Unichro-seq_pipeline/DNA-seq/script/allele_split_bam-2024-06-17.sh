#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
sample_barcode_file=info/2024-06-17.samples
locus_info=info/R10.locus_info2
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
ODIR=bowtie/20240617/$id

# Environment setting
ORGPATH=$( echo $PATH )
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1


cat $locus_info | awk '{print $1}' | while read snpid; do
   BAM=$ODIR/$snpid/$snpid.100bp_q30.bam
   
   # Get probe sequences for REF and ALT alleles
   ref_probe=$( cat info/rR10_snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="REF" ){print $3}}' )
   alt_probe=$( cat info/rR10_snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="ALT" ){print $3}}' )
   
   # Set environment
   export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH
   
   # Extract header from BAM file
   samtools view -H $BAM > $ODIR/$snpid/tmp.header
   
   # Process REF allele reads
   samtools view -F 64 $BAM |
     grep -F -v "XS:i:" - |
     awk -v ref_probe=$ref_probe 'BEGIN{FS="\t";OFS="\t"}{
       FLAG=64;
       $2=FLAG;
       SEQUENCE=$10;
       if(SEQUENCE ~ ref_probe){ print }
     }' |
     cat $ODIR/$snpid/tmp.header - |
     samtools view -bS - > $ODIR/$snpid/REF.bam

   samtools index $ODIR/$snpid/REF.bam
   
   # Process ALT allele reads
   samtools view -F 64 $BAM |
     grep -F -v "XS:i:" - |
     awk -v alt_probe=$alt_probe 'BEGIN{FS="\t";OFS="\t"}{
       FLAG=64;
       $2=FLAG;
       SEQUENCE=$10;
       if(SEQUENCE ~ alt_probe){ print }
     }' |
     cat $ODIR/$snpid/tmp.header - |
     samtools view -bS - > $ODIR/$snpid/ALT.bam

   samtools index $ODIR/$snpid/ALT.bam

   rm -f $ODIR/$snpid/tmp.header
done