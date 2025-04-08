#!/bin/sh
IFILE=$1
ODIR=$2
Nbase=$3
REF=${BOWTIE_REF:-"resources/GRCh38/GRCh38.primary_assembly"}
maxmismatch=2

echo "Written by K Ishigaki (2023-10-01)"
echo "Usage:"
echo "Parameter1: input file name containing two column data, an ID and a sequence"
echo "Parameter2: output directory name"
echo "Parameter3: this number of 3-end bases are used for mapping (or 'all' for full sequence)"

mkdir -p $ODIR

if [ "$Nbase" = "all" ]; then
    # If Nbase is "all", use the full sequence
    cat $IFILE | awk '{ print ">" $1; print $2 }' > $ODIR/input.fa
else
    # Otherwise, use only the last Nbase bases of the sequence
    cat $IFILE | awk -v Nbase=$Nbase '{
       print ">" $1;
       print substr($2, length($2) - Nbase + 1, Nbase)
    }' > $ODIR/input.fa
fi


bowtie --all -v $maxmismatch -x $REF -f $ODIR/input.fa > $ODIR/g38_map.org.txt

echo "ID Chr Position Strand N_mismatch Mismatch_info" |
perl -pe "s/ /\t/g" > $ODIR/g38_map.summary.txt

cat $ODIR/g38_map.org.txt |
awk 'BEGIN{FS="\t";OFS="\t"}{
   ID=$1;
   Strand=$2;
   Chr=$3;
   Position=$4;
   Mismatch_info=$8;
   N_mismatch=split(Mismatch_info,D,",");
   if( length(Mismatch_info)==0 ){
      N_mismatch=0;
      Mismatch_info="NA"
   };
   print ID, Chr, Position, Strand, N_mismatch, Mismatch_info
}' >> $ODIR/g38_map.summary.txt
