#!/bin/sh
jobid=${SGE_TASK_ID}

# Parameters
samplefile=info/2023-11-30.samples
indexinfo=info/tATAC.barcodes
id=$( cat -n $samplefile | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
ODIR=demultiplex/20231130/$id
F1=TMP/fastq/20231130/$id.preqc.umi.R1.fastq.gz
F2=TMP/fastq/20231130/$id.preqc.umi.R2.fastq.gz

# Environment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $ODIR

# Extract R1 reads with barcode
zcat $F1 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
awk 'BEGIN{FS="\t";OFS="\t"}{
   barcode = substr($2,18,8);
   print barcode, $1, $2, $3, $4
}' > $ODIR/tmp.R1

# Format R2 reads
zcat $F2 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}'  > $ODIR/tmp.R2

# Combine barcode, R1, and R2 information
paste $ODIR/tmp.R1 $ODIR/tmp.R2 > $ODIR/tmp.R1_R2

# Demultiplex reads based on barcode
cat $indexinfo | while read barcode; do
   # Process R1 reads matching the barcode
   cat $ODIR/tmp.R1_R2 |
   awk -v barcode=$barcode 'BEGIN{FS="\t";OFS="\t"}{
      if( $1 == barcode ){
         print $2 "\n" $3 "\n" $4 "\n" $5
      }
   }' |
   gzip -c - > $ODIR/$id.$barcode.preqc.umi.R1.fastq.gz
   
   # Process R2 reads matching the barcode
   cat $ODIR/tmp.R1_R2 |
   awk -v barcode=$barcode 'BEGIN{FS="\t";OFS="\t"}{
      if( $1 == barcode ){
         print $6 "\n" $7 "\n" $8 "\n" $9
      }
   }' |
   gzip -c - > $ODIR/$id.$barcode.preqc.umi.R2.fastq.gz
done  

# Create barcode-specific UMI files
cat $indexinfo | while read barcode; do
   F1=$ODIR/$id.$barcode.preqc.umi.R1.fastq.gz
   zcat $F1 |
   awk '{if(NR%4==2){print }}' |
   awk '{out = substr($1,1,17); print out}' |
   gzip -c > $ODIR/$id.$barcode.umi.gz
done