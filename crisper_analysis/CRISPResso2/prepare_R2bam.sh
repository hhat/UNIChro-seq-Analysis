#!/bin/bash
#

export PATH=/home/ha7477/tools/miniconda3/bin:${PATH}
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1
wd=/home/ha7477/works/umi/crispresso
eval "$(/home/ha7477/tools/miniconda3/bin/conda shell.bash hook)"
conda activate crispresso2

cd /home/imgkono/wd/img/crispr_qtl/DNAseq2
gRNA_file=/home/ha7477/works/umi/crispresso/data/R10_gRNA.txt

# Process each sample and SNP
{
  tail -n +2 /home/ha7477/works/umi/info/v1_DNAseq.txt
  tail -n +2 /home/ha7477/works/umi/info/v2_DNAseq.txt
} | while IFS=$'\t' read -r sample_raw replicate sample_id barcode bam_day; do
  tail -n +2 "${gRNA_file}" | while IFS=$'\t' read -r rezaid hg19pos hg38pos rsid jurkat guide strand; do
    snpid=$(echo "$hg38pos" | sed 's/^/chr/' | tr ':' '_')
    
    ref_probe=$( cat info/rR10_snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="REF" ){print $3}}' )
    alt_probe=$( cat info/rR10_snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="ALT" ){print $3}}' )
      
    old_bam="/home/imgkono/wd/img/crispr_qtl/DNAseq2/bowtie/${bam_day}/${sample_id}/${snpid}/${snpid}.100bp_q30_R2.bam"
    bam="/home/imgkono/wd/img/crispr_qtl/DNAseq2/bowtie/${bam_day}/${sample_id}/${snpid}/${snpid}.R2_after_probe.bam"
    
    samtools view -h $old_bam | \
    grep -F -v "XS:i:" | \
    awk -v ref_probe="$ref_probe" -v alt_probe="$alt_probe" '
      /^@/ {print; next}
      index($10, ref_probe) > 0 || index($10, alt_probe) > 0 {print}
    ' | samtools view -bS -o ${bam} -
    
    samtools index ${bam}
  done
done
