#!/bin/bash
#$ -S /bin/sh

wd=${1}
sample_id=${2}

fastq1=${wd}/trim_fastq/${sample_id}_val_1.fq.gz
fastq2=${wd}/trim_fastq/${sample_id}_val_2.fq.gz


outdir=${wd}/star/${sample_id}
mkdir -p $outdir

##Two pass mode(len 100)
/home/ha7477/tools/miniconda3/envs/de/bin/STAR \
--runMode alignReads \
--runThreadN 8 \
--genomeDir /home/ha7477/reference/rnaseq/star/grch38_gencodev26_GTExv8_99sjdb \
--twopassMode Basic \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFilterType BySJout \
--outFilterScoreMinOverLread 0.33 \
--outFilterMatchNmin 0 \
--outFilterMatchNminOverLread 0.33 \
--limitSjdbInsertNsj 1200000 \
--readFilesIn ${fastq1} ${fastq2} \
--readFilesCommand zcat \
--outFileNamePrefix ${outdir}/${sample_id}. \
--outSAMstrandField intronMotif \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--chimSegmentMin 15 \
--chimJunctionOverhangMin 15 \
--chimMainSegmentMultNmax 1 \
--chimOutJunctionFormat 0 \
--outSAMattributes NH HI AS nM NM ch \
--outSAMattrRGline ID:rg1 SM:sm1

