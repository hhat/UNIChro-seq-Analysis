#!/bin/bash
#$ -S /bin/sh

wd=${1}
sample_id=${2}

inputbam=${wd}/star/${sample_id}/${sample_id}.Aligned.sortedByCoord.out.bam

outdir=${wd}/rnaseqc2/${sample_id}
mkdir -p ${outdir}

genes_gtf=/home/ha7477/reference/rnaseq/GTExrepository/GENCODE_gencode.v26.GRCh38.genes.gtf
genome_fasta=/home/ha7477/reference/rnaseq/GTExrepository/references_Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta

##RNASeQC2
/home/ha7477/tools/tool/rnaseqc.v2.4.2.linux $genes_gtf ${inputbam} $outdir -s $sample_id -vv
