#!/bin/bash
#
# UMI Extraction and Processing Workflow
# Authors: Hiroaki Hatano / Kazuyoshi Ishigaki


wd=/home/ha7477/works/umi/stim_rs581
cd ${wd}


ORGPATH=$( echo $PATH )

#####################################################################
# 1. Experiment information files
#####################################################################

# Set data directory
dd=/home/imgkono/data/img/miseq/20231130/Fastq

# 1.1 Generate sample ID list
ls $dd |
grep "R1_001.fastq.gz" |
sed -e "s/_R1_001.fastq.gz//g" |
grep -v Undetermined > info/2023-11-30.samples

# 1.2 Create sample-barcode pair information
echo "01_S1_L001 CTCTCTAT
02_S2_L001 TATCCTCT
03_S3_L001 GTAAGGAG" > info/2023-11-30.sample_barcode_pair

#####################################################################
# 2. FASTQC
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 \
   ./script/fastqc-2023-11-30.sh

#####################################################################
# 3. Add UMI (extract first 17 bases from Read1)
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 \
   ./script/add_umi-2023-11-30.sh

#####################################################################
# 4. Demultiplex samples using custom barcodes
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 \
   ./script/demultiplex-2023-11-30.sh

#####################################################################
# 5. Adapter trimming with Cutadapt
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 \
   ./script/cutadapt-2023-11-30.sh

#####################################################################
# 6. Bowtie2 
#####################################################################

qsub -pe def_slot 4 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 \
   ./script/bowtie2-2023-11-30.sh


#####################################################################
# 7. Split BAM files by target regions
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 \
   ./script/split_bam_q30-2023-11-30.sh

#####################################################################
# 9. UMI counting
#####################################################################

batch=20231130
  qsub -pe def_slot 1 \
     -l s_vmem=64G \
     -cwd \
     -t 1:27 \
     ./script/umi_count_prep_R11_v1.sh ${batch}


#####################################################################
# 10. Filter UMI data 
#####################################################################

batch=20231130
  size_cutoff=300 #600
  count_cutoff=10
    qsub -pe def_slot 1 \
       -l s_vmem=4G \
       -cwd \
       ./script/filter_v12.sh ${size_cutoff} ${count_cutoff} ${batch}