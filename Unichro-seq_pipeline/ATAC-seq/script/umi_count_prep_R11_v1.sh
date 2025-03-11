#!/bin/sh
jobid=${SGE_TASK_ID}
batch=${1}

ORGPATH=$( echo $PATH )
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

locus_info=info/R11.locus_info2

batch2=`echo ${batch} | sed -E 's/([0-9]{4})([0-9]{2})([0-9]{2})/\1-\2-\3/'`
sample_barcode_file=info/${batch2}.sample_barcode_pair
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
fastqR1=/home/imgkono/data/img/novaseq/${batch}/Fastq/${id}_R1_001.fastq.gz
ODIR=bowtie/${batch}/${id}.${barcode}
demultiplex_dir=demultiplex/${batch}
Rscript=script/umi_count_v10.R

ALL_UMI_FILE="$ODIR/all_R1.umi.gz"
if [ ! -f "$ALL_UMI_FILE" ]; then
    echo "Creating $ALL_UMI_FILE"
    zcat $fastqR1 |
    awk '{if(NR%4==2){print }}' |
    awk '{out = substr($1,1,17); print out}' | gzip -c > "$ALL_UMI_FILE"
else
    echo "$ALL_UMI_FILE already exists, skipping creation"
fi

cat $locus_info | awk '{print $1}' | while read snpid; do
    BAM="$ODIR/$snpid/${snpid}.100bp_q30.bam"
    if [ -f "$BAM" ]; then
    echo "Processing $BAM"

    ref_probe=$( cat info/rR11_snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="REF" ){print $3}}' )
    alt_probe=$( cat info/rR11_snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="ALT" ){print $3}}' )

    export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

    # Process REF allele reads
    samtools view -f 146 $BAM |
    grep -F -v "XS:i:" - |
    awk 'BEGIN{FS="\t"}{
      Strand="-";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R1_1st_base=$8;
      size=$9;
      if(size< 0){size = -size};
      R2_1st_base = R1_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
    }' |
    grep -F $ref_probe > $ODIR/$snpid/$snpid.ref_umi_info

    samtools view -f 162 $BAM |
    grep -F -v "XS:i:" - |
    awk 'BEGIN{FS="\t"}{
      Strand="+";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R2_1st_base=$4;
      size=$9;
      if(size< 0){size = -size};
      R1_1st_base = R2_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
    }' |
    grep -F $ref_probe >> $ODIR/$snpid/$snpid.ref_umi_info

    # Process ALT allele reads
    samtools view -f 146 $BAM |
    grep -F -v "XS:i:" - |
    awk 'BEGIN{FS="\t"}{
      Strand="-";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R1_1st_base=$8;
      size=$9;
      if(size< 0){size = -size};
      R2_1st_base = R1_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
    }' |
    grep -F $alt_probe > $ODIR/$snpid/$snpid.alt_umi_info
   
    samtools view -f 162 $BAM |
    grep -F -v "XS:i:" - |
    awk 'BEGIN{FS="\t"}{
      Strand="+";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R2_1st_base=$4;
      size=$9;
      if(size< 0){size = -size};
      R1_1st_base = R2_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
    }' |
    grep -F $alt_probe >> $ODIR/$snpid/$snpid.alt_umi_info

    # Collect UMIs from other samples with same barcode
    find ${demultiplex_dir} -type f -name "*.${barcode}.umi.gz" ! -path "${demultiplex_dir}/${id}/*" -exec zcat {} + >> $ODIR/$snpid/${snpid}.samebarcode_R1.umi

    gzip $ODIR/$snpid/${snpid}.samebarcode_R1.umi

    # Run R script for UMI counting
    ORGPATH=$( echo $PATH )
    export PATH=/home/imgishi/miniconda3/envs/r4.2.2/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

    /home/ha7477/tools/miniconda3/envs/de/bin/R --vanilla --slave --args $ODIR/$snpid/$snpid.ref_umi_info $ODIR/$snpid/$snpid.alt_umi_info "$ALL_UMI_FILE" $ODIR/$snpid/${snpid}.samebarcode_R1.umi.gz < ${Rscript}

    else
        echo "Skipping $BAM - file not found"
    fi
done