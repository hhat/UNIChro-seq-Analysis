#!/bin/bash
#$ -S /bin/sh

chr=${1}

export PATH=/home/ha7477/tools/miniconda3/envs/atac1c/bin:/home/ha7477/tools/miniconda3/envs/atac1b/bin:/home/ha7477/tools/miniconda3/envs/atac1a/bin:${PATH}

cd /home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/b37


JAVA_TOOL_OPTIONS="-Xms48g -Xmx48g"

mkdir -p /home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/hg38liftover/rejected/

input=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/b37/common_ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EURcommon.vcf.gz
output=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/hg38liftover/common_ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_EURcommon_hg38liftover.vcf.gz


picard LiftoverVcf \
-Xms12g -Xmx12g \
MAX_RECORDS_IN_RAM=100000 \
WARN_ON_MISSING_CONTIG=true \
INPUT=${input} \
OUTPUT=${output} \
CHAIN=/home/ha7477/reference/chain/b37ToHg38.over.chain \
REJECT=/home/ha7477/share/to_imgkono2/reference/1KG_EURcommon/hg38liftover/rejected/${chr}_rejected_variants.vcf \
R=/home/ha7477/reference/rnaseq/GTExrepository/references_Homo_sapiens_assembly38_noALT_noHLA_noDecoy.fasta

