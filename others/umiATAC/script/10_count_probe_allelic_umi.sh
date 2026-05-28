#!/bin/bash
#$ -S /bin/bash

sample_id=${1}

wd=${wd:-/home/ha7477/works/umi/umiATAC}
probe_info=${probe_info:-/home/imgkono/wd/img/crispr_qtl/UNIChro_seq2/info/rTarget17_snp_probe_fw.info}
bam=${bam:-${wd}/umi_dedup_filtered/${sample_id}/${sample_id}_F780f2q30_autosomal_bfilt.bam}
odir=${odir:-${wd}/ASE_probe/${sample_id}}
snp_list=${snp_list:-"chr1_3900696_C_T chr1_16613699_G_A chr1_16904867_C_T chr1_201311493_G_A chr3_49340655_G_A chr4_118278629_G_C chr7_36152997_T_G chr7_91940815_C_T chr7_92621060_T_A chr7_139341719_G_T chr17_46194146_G_C chr17_57551724_G_A chr18_45687443_T_G"}

if [ -z "$sample_id" ]; then
  echo "Usage: $0 <sample_id>"
  exit 1
fi

ORGPATH=$( echo $PATH )
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:$ORGPATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p "$odir"

if [ ! -f "$bam" ]; then
  echo "BAM not found: $bam"
  exit 1
fi

for snpid in $snp_list; do
  ref_probe=$( awk -v snpid="$snpid" '{if( $1==snpid && $2=="REF" ){print $3}}' "$probe_info" )
  alt_probe=$( awk -v snpid="$snpid" '{if( $1==snpid && $2=="ALT" ){print $3}}' "$probe_info" )

  if [ -z "$ref_probe" ] || [ -z "$alt_probe" ]; then
    echo "Probe not found for $snpid, skipping"
    continue
  fi

  chr=$( echo "$snpid" | cut -d "_" -f 1 )
  pos=$( echo "$snpid" | cut -d "_" -f 2 )
  from=$( expr "$pos" - 100 )
  to=$( expr "$pos" + 100 )
  if [ "$from" -lt 1 ] ; then from=1 ; fi
  region="${chr}:${from}-${to}"

  mkdir -p "${odir}/${snpid}"
  echo "Processing $snpid (REF=$ref_probe, ALT=$alt_probe, region=$region)"

  samtools view -f 146 -F 256 "$bam" "$region" |
  grep -F -v "XS:i:" - |
  awk 'BEGIN{FS="\t"}{
    Strand="-";
    readid=$1;
    n=split(readid,D,"_");
    UMI=D[n];
    CIGAR=$6;
    gsub(/[0-9]/, "",CIGAR);
    R1_1st_base=$8;
    size=$9;
    if(size< 0){size = -size};
    R2_1st_base = R1_1st_base + size - 1;
    sequence=$10;
    print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
  }' |
  grep -F "$ref_probe" > "${odir}/${snpid}/${snpid}.ref_umi_info.tmp"

  samtools view -f 162 -F 256 "$bam" "$region" |
  grep -F -v "XS:i:" - |
  awk 'BEGIN{FS="\t"}{
    Strand="+";
    readid=$1;
    n=split(readid,D,"_");
    UMI=D[n];
    CIGAR=$6;
    gsub(/[0-9]/, "",CIGAR);
    R2_1st_base=$4;
    size=$9;
    if(size< 0){size = -size};
    R1_1st_base = R2_1st_base + size - 1;
    sequence=$10;
    print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
  }' |
  grep -F "$ref_probe" >> "${odir}/${snpid}/${snpid}.ref_umi_info.tmp"

  samtools view -f 146 -F 256 "$bam" "$region" |
  grep -F -v "XS:i:" - |
  awk 'BEGIN{FS="\t"}{
    Strand="-";
    readid=$1;
    n=split(readid,D,"_");
    UMI=D[n];
    CIGAR=$6;
    gsub(/[0-9]/, "",CIGAR);
    R1_1st_base=$8;
    size=$9;
    if(size< 0){size = -size};
    R2_1st_base = R1_1st_base + size - 1;
    sequence=$10;
    print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
  }' |
  grep -F "$alt_probe" > "${odir}/${snpid}/${snpid}.alt_umi_info.tmp"

  samtools view -f 162 -F 256 "$bam" "$region" |
  grep -F -v "XS:i:" - |
  awk 'BEGIN{FS="\t"}{
    Strand="+";
    readid=$1;
    n=split(readid,D,"_");
    UMI=D[n];
    CIGAR=$6;
    gsub(/[0-9]/, "",CIGAR);
    R2_1st_base=$4;
    size=$9;
    if(size< 0){size = -size};
    R1_1st_base = R2_1st_base + size - 1;
    sequence=$10;
    print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
  }' |
  grep -F "$alt_probe" >> "${odir}/${snpid}/${snpid}.alt_umi_info.tmp"

  awk -v pos="$pos" '{
    s = ($2 < $3) ? $2 : $3;
    e = ($2 > $3) ? $2 : $3;
    if(s <= pos && pos <= e) print
  }' "${odir}/${snpid}/${snpid}.ref_umi_info.tmp" > "${odir}/${snpid}/${snpid}.ref_umi_info"

  awk -v pos="$pos" '{
    s = ($2 < $3) ? $2 : $3;
    e = ($2 > $3) ? $2 : $3;
    if(s <= pos && pos <= e) print
  }' "${odir}/${snpid}/${snpid}.alt_umi_info.tmp" > "${odir}/${snpid}/${snpid}.alt_umi_info"

  rm -f "${odir}/${snpid}/${snpid}.ref_umi_info.tmp"
  rm -f "${odir}/${snpid}/${snpid}.alt_umi_info.tmp"

  ref_reads=$( wc -l < "${odir}/${snpid}/${snpid}.ref_umi_info" )
  alt_reads=$( wc -l < "${odir}/${snpid}/${snpid}.alt_umi_info" )
  ref_umi=$( awk '{print $1}' "${odir}/${snpid}/${snpid}.ref_umi_info" | sort -u | wc -l )
  alt_umi=$( awk '{print $1}' "${odir}/${snpid}/${snpid}.alt_umi_info" | sort -u | wc -l )

  echo -e "${sample_id}\t${snpid}\t${ref_reads}\t${alt_reads}\t${ref_umi}\t${alt_umi}" >> "${odir}/summary.txt"
  echo "  $snpid: REF_reads=$ref_reads ALT_reads=$alt_reads REF_umi=$ref_umi ALT_umi=$alt_umi"
done

echo "Done: ${sample_id}"
