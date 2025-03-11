#!/bin/bash
#$ -S /bin/sh

chr=${1}

export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}


export JAVA_HOME=/usr/local/package/java/current8
export PATH=/usr/local/package/java/current8/bin:$PATH
JAVA_TOOL_OPTIONS="-Xms50g -Xmx50g"

input=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/hg38liftover/common_ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EURcommon_hg38liftover.vcf.gz
output2=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/hg38liftover/common_ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EURcommon_hg38liftover_biallelicSNP.vcf.gz
output3=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/hg38liftover/common_ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EURcommon_hg38liftover_biallelicSNP_uniq.vcf

/home/ha7477/tools/tool/gatk-4.2.0.0/gatk SelectVariants -R /home/ha7477/reference/rnaseq/GTExrepository/references_Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta -V ${input} -O ${output2} --restrict-alleles-to BIALLELIC --select-type-to-include SNP

zcat ${output2} | awk '!colname[$1,$2]++{print}' > ${output3}
bgzip ${output3}
tabix -p vcf ${output3}.gz
