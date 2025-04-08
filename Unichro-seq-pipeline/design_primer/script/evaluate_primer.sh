#!/bin/sh
ODIR=$1
script_dir=$(dirname "$0")

 #prepare input primer data
cat $ODIR/target.primers |
awk 'BEGIN{FS="\t"}{ if(NR!=1){
   print $1, $4
}}' > $ODIR/tmp.input

 #mapping using all primer sequences
${script_dir}/short_read_map_g38.sh $ODIR/tmp.input $ODIR "all"
mv $ODIR/g38_map.org.txt $ODIR/g38_map.org.allseq.txt
mv $ODIR/g38_map.summary.txt $ODIR/g38_map.allseq.summary.txt

 #mapping using 3' end 15bp of primer sequences
${script_dir}/short_read_map_g38.sh $ODIR/tmp.input $ODIR 15
mv $ODIR/g38_map.org.txt $ODIR/g38_map.org.15.txt
mv $ODIR/g38_map.summary.txt $ODIR/g38_map.15.summary.txt

rm -f $ODIR/tmp.input
