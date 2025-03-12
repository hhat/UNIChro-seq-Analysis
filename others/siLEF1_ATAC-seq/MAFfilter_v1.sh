#!/bin/bash
#$ -S /bin/sh

chr=${1}

export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}


input=/home/ha7477/reference/1KG/vcf/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
output=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/b37/common_ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EURcommon.vcf
bcftools norm -m - ${input} | bcftools view -i "EUR_AF >0.05 && EUR_AF < 0.95" > ${output}
bgzip ${output}
tabix -p vcf ${output}.gz
