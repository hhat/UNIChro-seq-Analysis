#!/bin/bash
#
# DNAseq Analysis Workflow
# Authors: Michihiro Kono / Kazuyoshi Ishigaki / Hiroaki Hatano




wd=/home/imgkono/wd/img/crispr_qtl/DNAseq2/
dd=/home/imgkono/data/img/miseq/20240617/Fastq
cd $wd


ORGPATH=$( echo $PATH )
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH

#####################################################################
# 1. Experiment information files
#####################################################################

# Generate sample ID list
ls $dd |
grep "R1_001.fastq.gz" |
sed -e "s/_R1_001.fastq.gz//g" |
grep -v Undetermined > info/2024-06-17.samples

#####################################################################
# 2. FASTQC
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:16 \
   ./script/fastqc-2024-06-17.sh

#####################################################################
# 3. Adapter trimming with Cutadapt
#####################################################################

# Submit cutadapt job array
qsub -pe def_slot 1 \
   -l s_vmem=20G,mem_req=20G \
   -cwd \
   -t 1:16 \
   ./script/cutadapt-2024-06-17.sh

# Generate read count summary
mkdir -p analysis/2024-06-17/read_count

OFILE=analysis/2024-06-17/read_count/cutadapt_round1_readcount.txt
echo "ID N_Input_R1 N_Input_R2 N_Output_R1 N_Output_R2" > $OFILE
cat info/2024-06-17.samples | while read line; do
   id=$( echo $line | awk '{print $1}' )
   ODIR=cutadapt/20240617/$id
   cat $ODIR/round1.count >> $OFILE
done


#####################################################################
# 5. Bowtie2
#####################################################################

qsub -pe def_slot 4 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:16 \
   ./script/bowtie2-2024-06-17.sh

#####################################################################
# 6. Split BAM files by target regions
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:16 \
   ./script/split_bam_q30-2024-06-17.sh

#####################################################################
# 7. Allele-specific read extraction
#####################################################################

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:16 \
   ./script/allele_split_bam-2024-06-17.sh

#####################################################################
# 8. Allele counting from split BAM files
#####################################################################

# Count reads for each allele
mkdir -p analysis/2024-06-17/allelic_ratio
OFILE=analysis/2024-06-17/allelic_ratio/alelle_count_bam.txt
echo "ID SNP REF_count ALT_count" > $OFILE

cat info/2024-06-17.samples | while read id; do
  cat info/R10.locus_info2 | awk '{print $1}' | while read snpid; do
    echo "Processing $id $snpid"
    
    REF_BAM=bowtie/20240617/$id/$snpid/REF.bam
    ALT_BAM=bowtie/20240617/$id/$snpid/ALT.bam
    REF_count=$( samtools view -c $REF_BAM )
    ALT_count=$( samtools view -c $ALT_BAM )
    echo $id $snpid $REF_count $ALT_count >> $OFILE
  done
done
