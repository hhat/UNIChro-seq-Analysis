#!/bin/bash
#$ -S /bin/sh

sample_id=${1}
chr=${2}

export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

wd=/home/ha7477/share/to_imgkono2/result/Jurkat
ref=/home/ha7477/reference/rnaseq/GTExrepository/references_Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta 


dd=${wd}/ASE/${sample_id}
mkdir -p ${dd}

bam=${wd}/bowtie2/${sample_id}/chr/chr${chr}_${sample_id}_markdup_F1804f2q30_autosomal_bfilt_addRG.bam
region=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/hg38liftover/common_ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EURcommon_hg38liftover_biallelicSNP_uniq.vcf.gz

#bcftools mpileup -a "AD,DP" -Ov -f ${ref} ${bam} --max-depth 10000 -R ${region} > ${dd}/${sample_id}_step1.vcf
bcftools mpileup -a "AD,DP,DP4" -Ov -f ${ref} ${bam} --max-depth 10000 -R ${region} > ${dd}/chr${chr}_${sample_id}_step1.vcf


bgzip ${dd}/chr${chr}_${sample_id}_step1.vcf
tabix -f -p vcf ${dd}/chr${chr}_${sample_id}_step1.vcf.gz

dpcut=0
echo $sample_id > ${dd}/${sample_id}_chr${chr}.name
bcftools call -mA -Ou ${dd}/chr${chr}_${sample_id}_step1.vcf.gz | bcftools view -i "FORMAT/DP>${dpcut}" | bcftools norm -m -  | bcftools reheader -s ${dd}/${sample_id}_chr${chr}.name > ${dd}/chr${chr}_${sample_id}_step2.vcf
rm  ${dd}/${sample_id}_chr${chr}.name

bgzip ${dd}/chr${chr}_${sample_id}_step2.vcf
tabix -f -p vcf ${dd}/chr${chr}_${sample_id}_step2.vcf.gz
