#!/bin/bash
#$ -S /bin/sh


export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

cpu=8

sample_id=${1}

saf=/home/ha7477/share/to_imgkono2/result/siLEF1_ATAC/merged_peak/siLEF1_ATAC_all_peaks.narrowPeak_homerannotation.bed.saf

tarbam=/home/ha7477/share/to_imgkono2/result/siLEF1_ATAC/bowtie2/${sample_id}/${sample_id}_markdup_F1804f2q30_autosomal_bfilt.bam

featureCounts -F SAF -T $cpu -p -a $saf -o /home/ha7477/share/to_imgkono2/result/siLEF1_ATAC/featurecounts/${sample_id}.featureCounts.txt ${tarbam}

