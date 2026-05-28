#!/bin/bash
#$ -S /bin/sh


export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

wd=${1}
sample_id=${2}
core=8



##Mark Duplicate 
tmpfullbam=${wd}/bowtie2/${sample_id}/${sample_id}_sorted.bam


JAVA_TOOL_OPTIONS="-Xms12g -Xmx12g"
picard MarkDuplicates I=${tmpfullbam} O=${tmpfullbam}_markdup.bam M=${tmpfullbam}_metrics.txt


samtools flagstat $tmpfullbam > ${tmpfullbam}_flagstat

##samtools 
samtools view -F 1804 -f 2 -q 30 ${tmpfullbam}_markdup.bam -Obam | samtools sort - -o ${tmpfullbam}_markdup_F1804f2q30.bam

rm ${tmpfullbam}_markdup.bam

inputbam=`basename ${tmpfullbam}_markdup_F1804f2q30.bam`
sd=`dirname ${tmpfullbam}_markdup_F1804f2q30.bam`
cd $sd

samtools flagstat $inputbam > ${inputbam}_flagstat

BLACKLIST=/home/ha7477/reference/atac/ATACblacklist/hg38.blacklist.bed.gz
outputbam=${sample_id}_markdup_F1804f2q30_autosomal_bfilt.bam

samtools index ${inputbam}
samtools view -o - ${inputbam} `seq 1 22 | sed 's/^/chr/'` -b | bedtools intersect -nonamecheck -v -abam "stdin" -b ${BLACKLIST} | samtools sort -@ ${core} > ${outputbam}
samtools index ${outputbam}

samtools flagstat $outputbam > ${outputbam}_flagstat

macs2 callpeak -t ${outputbam} -g hs -f BAMPE -n ${sample_id} -B -q 0.05

TSS=/home/ha7477/reference/atac/TSS/refseq_hg38_UCSC_TSS.bed

ataqv \
--name ${sample_id} --tss-file $TSS \
--less-redundant \
--peak-file ${sample_id}_peaks.narrowPeak \
--metrics-file ${sample_id}.ataqv.json --ignore-read-groups human ${outputbam} > ${sample_id}.ataqv.log

mkarv ${sample_id} ${sample_id}.ataqv.json

rm ${sample_id}_control_lambda.bdg
rm ${sample_id}_treat_pileup.bdg
rm ${sample_id}_summits.bed
rm ${inputbam}
rm ${inputbam}.bai
