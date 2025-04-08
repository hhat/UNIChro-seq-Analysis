#!/bin/sh
IFILE=$1
ODIR=$2
FASTA=${GENOME_FASTA:-"resources/GRCh38.primary_assembly.genome.fa.gz"}

echo "Written by K Ishigaki (2023-10-01)"
echo "Usage:"
echo "Parameter1: input file name containing three column data: snpid, region, fw_1_rv_2"
echo "Parameter2: output directory name"

 #1, Fw sequence
cat $IFILE |
awk 'BEGIN{OFS="\t"}{
   if($3 == 1){
      split($1, D, "_");
      print D[1] ":" D[2] - 199 "-" D[2] , $1
   }
}' |
sort -k 1b,1 | uniq > $ODIR/tmp.target.region.info

cut -f 1 $ODIR/tmp.target.region.info |
samtools faidx $FASTA --length 200 --region-file - |
awk '{
   if( NR % 2 == 1 ){ 
      gsub(">","",$1); 
      printf $1 "\t"
   } else {
      print
   }
}' |
sort -k 1b,1 |
join $ODIR/tmp.target.region.info - |
awk '{print ">" $2; print $3}' > $ODIR/tmp.target.fw.fa
    #the last base is the target snp position

 #2, Rv sequence
cat $IFILE |
awk 'BEGIN{OFS="\t"}{
   if($3 == 2){
      split($1, D, "_");
      print D[1] ":" D[2] "-" D[2] + 199, $1
   }
}' |
sort -k 1b,1 | uniq > $ODIR/tmp.target.region.info

cut -f 1 $ODIR/tmp.target.region.info |
samtools faidx $FASTA --reverse-complement --length 200 --region-file - |
awk '{
   if( NR % 2 == 1 ){ 
      gsub(">","",$1); 
      gsub("/rc","",$1); 
      printf $1 "\t"
   } else {
      print
   }
}' |
sort -k 1b,1 |
join $ODIR/tmp.target.region.info - |
awk '{print ">" $2; print $3}' > $ODIR/tmp.target.rv.fa
    #the last base is the target snp position (RC)

 #3, merge data
cat $ODIR/tmp.target.fw.fa $ODIR/tmp.target.rv.fa > $ODIR/target.200bp.fa

rm -f $ODIR/tmp*
