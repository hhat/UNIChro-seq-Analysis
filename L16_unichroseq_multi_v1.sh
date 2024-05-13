 #Shirokane-IMG
 #2023-12-30
 #Project-K
 #Sample1-3 Jurkat 17 target
 #/home/imgkono/data/img/miseq/20231130/Fastq
 #Modified point
 #1, prepared scripts applicable to many projects
 #2, custom UMI counting

ssh -L 1234:localhost:1234 shirogane2

qlogin -l os7 -l s_vmem=20G
unset PROMPT_COMMAND

wd=/home/imgishi/wd/crisprqtl
cd $wd
mkdir -p log
mkdir -p jupyter

 #eval "$(/home/imgishi/miniconda3/bin/conda shell.bash hook)"
 #conda activate genetics1
 # /home/imgishi/miniconda3/envs/jupyter/bin/jupyter lab --port 1234 --no-browser

 #save the current path
ORGPATH=$( echo $PATH )





#####0, Experiment information files

dd=/home/imgkono/data/img/miseq/20231130/Fastq

 #1, sample id
ls $dd |
grep "R1_001.fastq.gz" |
sed -e "s/_R1_001.fastq.gz//g" |
grep -v Undetermined > info/2023-11-30.samples

 #check file
cat info/2023-11-30.samples #ok

 #2, barcode list: already exist
 #echo "CTCTCTAT
 #TATCCTCT
 #GTAAGGAG
 #ACTGCATA
 #AAGGAGTA
 #CTAAGCCT
 #CGTCTAAT
 #TCTCTCCG
 #TCGACTAG
 #TTCTAGCT
 #CCTAGAGT
 #GCGTAAGA
 #CTATTAAG
 #AAGGCTAT
 #GAGCCTTA
 #TTATGCGA" > info/tATAC.barcodes

 #3, correct sample id - barcode pairs
echo "01_S1_L001 CTCTCTAT
02_S2_L001 TATCCTCT
03_S3_L001 GTAAGGAG" > info/2023-11-30.sample_barcode_pair






#####1, FASTQC (pre-QC)

qsub  -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/fastqc-2023-11-30.log \
   -e log/fastqc-2023-11-30.error \
   ./script/fastqc-2023-11-30.sh
    #Your job-array 84357642.1-3:1 ("fastqc-2023-11-30.sh") has been submitted

 #check output: ok
cat info/2023-11-30.samples | while read id; do
   ll fastqc/20231130/$id.R1/${id}_R1_001_fastqc.html
   ll fastqc/20231130/$id.R2/${id}_R2_001_fastqc.html
done

 # ./script/fastqc-2023-11-30.sh
 #----------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
samplefile=info/2023-11-30.samples
mainodir=fastqc/20231130
dd=/home/imgkono/data/img/miseq/20231130/Fastq
id=$( cat -n $samplefile | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
F1=$dd/${id}_R1_001.fastq.gz
F2=$dd/${id}_R2_001.fastq.gz
 #####################

 #Enviroment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

 #R1
odir=$mainodir/$id.R1
mkdir -p $odir
fastqc -o $odir -f fastq $F1

 #R2
odir=$mainodir/$id.R2
mkdir -p $odir
fastqc -o $odir -f fastq $F2

 #----------------------------------------#






#####2, Add UMI (the first 17 bases)

qsub  -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/add_umi-2023-11-30.log \
   -e log/add_umi-2023-11-30.error \
   ./script/add_umi-2023-11-30.sh
    #Your job-array 84357692.1-3:1 ("add_umi-2023-11-30.sh") has been submitted

 #check output: ok!
cat info/2023-11-30.samples | while read id; do
   ll TMP/fastq/20231130/$id.preqc.umi.R1.fastq.gz
   ll TMP/fastq/20231130/$id.preqc.umi.R2.fastq.gz
   ll TMP/fastq/20231130/$id.umi.gz
done

 # ./script/add_umi-2023-11-30.sh
 #----------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
samplefile=info/2023-11-30.samples
mainodir=TMP/fastq/20231130
dd=/home/imgkono/data/img/miseq/20231130/Fastq
id=$( cat -n $samplefile | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
F1=$dd/${id}_R1_001.fastq.gz
F2=$dd/${id}_R2_001.fastq.gz
 #####################

 #Enviroment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $mainodir

 #get umi (first 17 bases in Read1)
zcat $F1 |
awk '{if(NR%4==2){print }}' |
awk '{out = substr($1,1,17); print out}' > $mainodir/$id.umi

 #add UMP to read name (R1)
zcat $F1 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
paste $mainodir/$id.umi - |
awk 'BEGIN{FS="\t"}{
   umi=$1;
   readid=$2;
   sequence=$3;
   strand=$4;
   quality=$5;
   split(readid, D, " ");
   print D[1] "_" umi " " D[2] "\n" sequence "\n" strand "\n" quality
}' |
gzip -c - > $mainodir/$id.preqc.umi.R1.fastq.gz

 #add UMI to read name (R2)
zcat $F2 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
paste $mainodir/$id.umi - |
awk 'BEGIN{FS="\t"}{
   umi=$1;
   readid=$2;
   sequence=$3;
   strand=$4;
   quality=$5;
   split(readid, D, " ");
   print D[1] "_" umi " " D[2] "\n" sequence "\n" strand "\n" quality
}' |
gzip -c - > $mainodir/$id.preqc.umi.R2.fastq.gz

gzip -f $mainodir/$id.umi

 #----------------------------------------#







#####3, Demultiplex
 #using custom barcode
 # requiring perfect match

qsub  -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/demultiplex-2023-11-30.log \
   -e log/demultiplex-2023-11-30.error \
   ./script/demultiplex-2023-11-30.sh
    #Your job-array 84357816.1-3:1 ("demultiplex-2023-11-30.sh") has been submitted

 #check output: ok!
cat info/2023-11-30.samples | while read id; do
cat info/tATAC.barcodes | while read barcode;do
   ll demultiplex/20231130/$id/$id.$barcode.preqc.umi.R1.fastq.gz
   ll demultiplex/20231130/$id/$id.$barcode.preqc.umi.R2.fastq.gz
   ll demultiplex/20231130/$id/$id.$barcode.umi.gz
done
done

 #Read count: ok
mkdir -p analysis/2023-11-30/read_count
OFILE=analysis/2023-11-30/read_count/preqc_readcount_demultiplex.txt
echo "ID Barcode N_Input_R1 N_Input_R2 N_Output_R1 N_Output_R2" > $OFILE

cat info/2023-11-30.samples | while read id;do
cat info/tATAC.barcodes | while read barcode;do
   echo $id : $barcode ...
   
   F1=TMP/fastq/20231130/$id.preqc.umi.R1.fastq.gz
   F2=TMP/fastq/20231130/$id.preqc.umi.R2.fastq.gz
   OF1=demultiplex/20231130/$id/$id.$barcode.preqc.umi.R1.fastq.gz
   OF2=demultiplex/20231130/$id/$id.$barcode.preqc.umi.R2.fastq.gz
   
   N_Input_R1=$( zcat $F1 | wc -l | awk '{ print $1 /4 }' )
   N_Input_R2=$( zcat $F2 | wc -l | awk '{ print $1 /4 }' )
   N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
   N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )
   
   echo $id $barcode $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 >> $OFILE
   
done
done

 #check: ok
cat analysis/2023-11-30/read_count/preqc_readcount_demultiplex.txt |
awk '{if($5 > 10000){print}}'

 # ./script/demultiplex-2023-11-30.sh
 #----------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
samplefile=info/2023-11-30.samples
indexinfo=info/tATAC.barcodes
id=$( cat -n $samplefile | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
ODIR=demultiplex/20231130/$id
F1=TMP/fastq/20231130/$id.preqc.umi.R1.fastq.gz
F2=TMP/fastq/20231130/$id.preqc.umi.R2.fastq.gz
 #####################

 #Enviroment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $ODIR

 #R1 with cell barcode
zcat $F1 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}' |
awk 'BEGIN{FS="\t";OFS="\t"}{
   barcode = substr($2,18,8);
   print barcode, $1, $2, $3, $4
}' > $ODIR/tmp.R1

 #R2
zcat $F2 |
awk '{ 
   printf("%s",$0);
   if(NR%4==0){printf("\n")} else { printf("\t")}
}'  > $ODIR/tmp.R2

 #barcode, R1, R2
paste $ODIR/tmp.R1 $ODIR/tmp.R2 > $ODIR/tmp.R1_R2

 #demultiplex
cat $indexinfo | while read barcode;do
   
   echo Starting demultiplexing by $barcode from $id sample ...
   
   #R1
   cat $ODIR/tmp.R1_R2 |
   awk -v barcode=$barcode 'BEGIN{FS="\t";OFS="\t"}{
      if( $1 == barcode ){
         print $2 "\n" $3 "\n" $4 "\n" $5
      }
   }' |
   gzip -c - > $ODIR/$id.$barcode.preqc.umi.R1.fastq.gz
   
   #R2
   cat $ODIR/tmp.R1_R2 |
   awk -v barcode=$barcode 'BEGIN{FS="\t";OFS="\t"}{
      if( $1 == barcode ){
         print $6 "\n" $7 "\n" $8 "\n" $9
      }
   }' |
   gzip -c - > $ODIR/$id.$barcode.preqc.umi.R2.fastq.gz
   
done  

 # Barcode-specific TMP_UMI file
cat $indexinfo | while read barcode;do
   F1=$ODIR/$id.$barcode.preqc.umi.R1.fastq.gz
   zcat $F1 |
   awk '{if(NR%4==2){print }}' |
   awk '{out = substr($1,1,17); print out}' |
   gzip -c > $ODIR/$id.$barcode.umi.gz
   
done

 #----------------------------------------#








#####4, Cutadapt

qsub  -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/cutadapt-2023-11-30.log \
   -e log/cutadapt-2023-11-30.error \
   ./script/cutadapt-2023-11-30.sh
    #Your job-array 84358150.1-3:1 ("cutadapt-2023-11-30.sh") has been submitted

 #check outputs: ok
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   ll cutadapt/20231130/$id/$id.$barcode.postqc.umi.R1.fastq.gz
   ll cutadapt/20231130/$id/$id.$barcode.postqc.umi.R2.fastq.gz
done

 #Read count summary
mkdir -p analysis/2023-11-30/read_count

OFILE=analysis/2023-11-30/read_count/cutadapt_round1_readcount.txt
echo "ID Barcode N_Input_R1 N_Input_R2 N_Output_R1 N_Output_R2 N_Output_notrim_R1 N_Output_notrim_R2" > $OFILE
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   ODIR=cutadapt/20231130/$id
   cat $ODIR/round1.count >> $OFILE
done

OFILE=analysis/2023-11-30/read_count/cutadapt_round2_readcount.txt
echo "ID Barcode N_Input_R1 N_Input_R2 N_Output_R1 N_Output_R2 N_Output_notrim_R1 N_Output_notrim_R2" > $OFILE
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   ODIR=cutadapt/20231130/$id
   cat $ODIR/round2.count >> $OFILE
done

OFILE=analysis/2023-11-30/read_count/cutadapt_round3_readcount.txt
echo "ID Barcode N_Input_R1 N_Input_R2 N_Output_R1 N_Output_R2" > $OFILE
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   ODIR=cutadapt/20231130/$id
   cat $ODIR/round3.count >> $OFILE
done

 # ./script/cutadapt-2023-11-30.sh
 #----------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
sample_barcode_file=info/2023-11-30.sample_barcode_pair
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
IDIR=demultiplex/20231130/$id
ODIR=cutadapt/20231130/$id
 #####################

 #Enviroment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p $ODIR


 ### Round1 ###
 #Trim AGATGTGTATAAGAGACAG from 5' ONLY ONCE
 #Outputs reads without AGATGTGTATAAGAGACAG to diagnosis the experiments
IF1=$IDIR/$id.$barcode.preqc.umi.R1.fastq.gz
IF2=$IDIR/$id.$barcode.preqc.umi.R2.fastq.gz
OF1=$ODIR/$id.$barcode.postqc_tmpA.umi.R1.fastq.gz
OF2=$ODIR/$id.$barcode.postqc_tmpA.umi.R2.fastq.gz
OF1_notrim=$ODIR/$id.$barcode.postqc_tmpA_notrim.umi.R1.fastq.gz
OF2_notrim=$ODIR/$id.$barcode.postqc_tmpA_notrim.umi.R2.fastq.gz

cutadapt \
   -g "AGATGTGTATAAGAGACAG;min_overlap=19" \
   --times 1 \
   --untrimmed-output  $OF1_notrim \
   --untrimmed-paired-output  $OF2_notrim \
   -o  $OF1 -p  $OF2 \
   $IF1 $IF2

N_Input_R1=$( zcat $IF1 | wc -l | awk '{ print $1 /4 }' )
N_Input_R2=$( zcat $IF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R1=$( zcat $OF1_notrim | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R2=$( zcat $OF2_notrim | wc -l | awk '{ print $1 /4 }' )

echo $id $barcode $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 $N_Output_notrim_R1 $N_Output_notrim_R2 >> $ODIR/round1.count


 ### Round2 ###
 #Trim AGATGTGTATAAGAGACAG from 5'MULTIPLE TIMES
 #In this round, the reads WITHOUT triming is the ones we should use in the downstream analyses.
IF1=$ODIR/$id.$barcode.postqc_tmpA.umi.R1.fastq.gz
IF2=$ODIR/$id.$barcode.postqc_tmpA.umi.R2.fastq.gz
OF1=$ODIR/$id.$barcode.postqc_tmpB.umi.R1.fastq.gz
OF2=$ODIR/$id.$barcode.postqc_tmpB.umi.R2.fastq.gz
OF1_notrim=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R1.fastq.gz
OF2_notrim=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R2.fastq.gz

cutadapt \
   -g "AGATGTGTATAAGAGACAG;min_overlap=19" \
   --times 5 \
   --untrimmed-output  $OF1_notrim \
   --untrimmed-paired-output  $OF2_notrim \
   -o  $OF1 -p  $OF2 \
   $IF1 $IF2

N_Input_R1=$( zcat $IF1 | wc -l | awk '{ print $1 /4 }' )
N_Input_R2=$( zcat $IF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R1=$( zcat $OF1_notrim | wc -l | awk '{ print $1 /4 }' )
N_Output_notrim_R2=$( zcat $OF2_notrim | wc -l | awk '{ print $1 /4 }' )

echo $id $barcode $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 $N_Output_notrim_R1 $N_Output_notrim_R2 >> $ODIR/round2.count


 ### Round3 ###
 #Trim CTGTCTCTTATACACATCT from 3' MULTIPLE TIMES
IF1=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R1.fastq.gz
IF2=$ODIR/$id.$barcode.postqc_tmpB_notrim.umi.R2.fastq.gz
OF1=$ODIR/$id.$barcode.postqc.umi.R1.fastq.gz
OF2=$ODIR/$id.$barcode.postqc.umi.R2.fastq.gz

cutadapt \
   -a CTGTCTCTTATACACATCT \
   -A CTGTCTCTTATACACATCT \
   --minimum-length 20 \
   --times 5 \
   -o  $OF1 -p  $OF2 \
   $IF1 $IF2

N_Input_R1=$( zcat $IF1 | wc -l | awk '{ print $1 /4 }' )
N_Input_R2=$( zcat $IF2 | wc -l | awk '{ print $1 /4 }' )
N_Output_R1=$( zcat $OF1 | wc -l | awk '{ print $1 /4 }' )
N_Output_R2=$( zcat $OF2 | wc -l | awk '{ print $1 /4 }' )

echo $id $barcode $N_Input_R1 $N_Input_R2 $N_Output_R1 $N_Output_R2 > $ODIR/round3.count

 #MEMO:
 #AGATGTGTATAAGAGACAG (R1,5'): -g (ME in a positive strand; in a standard library, this is not sequence because this works as a sequence primer, but in our pipeline, we do not use illumina primer cocktail)
 #CTGTCTCTTATACACATCT (R1,3'): -a (ME in a revers strand)
 #CTGTCTCTTATACACATCT (R2,3'): -A (ME in a revers strand)

 #----------------------------------------#





#####5, FASTQC (post-QC)

qsub  -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/fastqc_postqc-2023-11-30.log \
   -e log/fastqc_postqc-2023-11-30.error \
   ./script/fastqc_postqc-2023-11-30.sh

 #check output: ok
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   ll fastqc/20231130/$id.R1.postqc/$id.$barcode.postqc.umi.R1_fastqc.zip
   ll fastqc/20231130/$id.R2.postqc/$id.$barcode.postqc.umi.R2_fastqc.zip
done

 # ./script/fastqc_postqc-2023-11-30.sh
 #----------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
mainodir=fastqc/20231130
sample_barcode_file=info/2023-11-30.sample_barcode_pair
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
ODIR=cutadapt/20231130/$id
F1=$ODIR/$id.$barcode.postqc.umi.R1.fastq.gz
F2=$ODIR/$id.$barcode.postqc.umi.R2.fastq.gz
 #####################

 #Enviroment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

 #R1
odir=$mainodir/$id.R1.postqc
mkdir -p $odir
fastqc -o $odir -f fastq $F1

 #R2
odir=$mainodir/$id.R2.postqc
mkdir -p $odir
fastqc -o $odir -f fastq $F2

 #----------------------------------------#




#####6, Bowtie2

 #4 threads job
qsub -pe def_slot 4 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/bowtie2-2023-11-30.log \
   -e log/bowtie2-2023-11-30.error \
   ./script/bowtie2-2023-11-30.sh
     #Your job-array 84414937.1-3:1 ("bowtie2-2023-11-30.sh") has been submitted
     #MEMO:
     #-q INT   only include reads with mapping quality >= INT [0]
     #bowtie2 mapping quality: higher = more unique
     #Aligners characterize their degree of confidence in the point of origin by reporting a mapping quality: a non-negative integer Q = -10 log10 p, where p is an estimate of the probability that the alignment does not correspond to the read's true point of origin. Mapping quality is sometimes abbreviated MAPQ, and is recorded in the SAM MAPQ field.

 #ouput check: ok
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   ll bowtie/20231130/$id.$barcode/d1.sorted.bam
done

 #Add q 30 filtering (updated)
 #many indels maybe due to low quality
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   IBAM=bowtie/20231130/$id.$barcode/d1.sorted.bam
   OBAM=bowtie/20231130/$id.$barcode/d1.sorted.q30.bam
   
   samtools view -f 2 -q 30 $IBAM -o $OBAM
   
   samtools index $OBAM
   
done



 # ./script/bowtie2-2023-11-30.sh
 #-------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
threads=4 #number of alignment threads to launch
sample_barcode_file=info/2023-11-30.sample_barcode_pair
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
F1=cutadapt/20231130/$id/$id.$barcode.postqc.umi.R1.fastq.gz
F2=cutadapt/20231130/$id/$id.$barcode.postqc.umi.R2.fastq.gz
ODIR=bowtie/20231130/$id.$barcode
OFILE=bowtie/20231130/$id.$barcode/d1
Index=/home/imgishi/reference/Bowtie2_ref/GRCh38/GRCh38
 #####################

 #Enviroment setting
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

 #Mapping
mkdir -p $ODIR

bowtie2 -x $Index \
   -1 $F1 \
   -2 $F2 \
   --threads $threads \
   --no-discordant \
   --maxins  1000 \
   -S $OFILE.sam
    #comment: added --maxins  1000

 #Filter mapping outputs
samtools view -bS $OFILE.sam > $OFILE.bam

samtools view -f 2 -q 10 $OFILE.bam -o $OFILE.qc.bam
    #-f 2: read mapped in proper pair

samtools sort $OFILE.qc.bam -o $OFILE.sorted.bam

samtools index $OFILE.sorted.bam

rm -f $OFILE.sam
rm -f $OFILE.qc.bam

 #-------------------------------------# 





#####7-1, Locus inforamtion (updated)

 #input file
 #col1: snp id (GRCh38; chr_pos_ref_alt)
 #col2: inner primer sequence
 #col3: prmer on forward (1) or reverse strand (2)
 #./info/2023-11-30.locus_info1

head -n 3 ./info/2023-11-30.locus_info1
    #chr1_16613699_G_A   ACTCGGTCCTTCCCCTGGGTTACTC  2
    #chr1_16904867_C_T   CCCCACCCTCCTCCACTCAG              2
    #chr1_201311493_G_A TTTGAGCTCCCAGAAAGACTGGGTG  2

export PATH=/home/imgishi/miniconda3/envs/r4.2.2/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

R
library(Biostrings)

d1 <- read.table("info/2023-11-30.locus_info1",header=F)
head(d1,n=3)
    #chr1_16613699_G_A      ACTCGGTCCTTCCCCTGGGTTACTC    2
    #chr1_16904867_C_T      CCCCACCCTCCTCCACTCAG                2
    #chr1_201311493_G_A     TTTGAGCTCCCAGAAAGACTGGGTG  2

d1$RC <- sapply(as.character(d1$V2),function(x){
    return(as.character(reverseComplement(DNAString(x))))
}) #add reverse complement sequence to the 4th column

OFILE="info/2023-11-30.locus_info2"
write.table(d1,OFILE,col.names=F,row.names=F,append=F,quote=F)
  #R off

head -n 3 ./info/2023-11-30.locus_info2
   #e.g.,
   #chr1_16613699_G_A    ACTCGGTCCTTCCCCTGGGTTACTC   2   GAGTAACCCAGGGGAAGGACCGAGT
   #chr1_16904867_C_T     CCCCACCCTCCTCCACTCAG              2   CTGAGTGGAGGAGGGTGGGG
   #chr1_201311493_G_A   TTTGAGCTCCCAGAAAGACTGGGTG  2   CACCCAGTCTTTCTGGGAGCTCAAA

   #col1: snp id (GRCh38; chr_pos_ref_alt)
   #col2: inner primer sequence
   #col3: prmer on forward (1) or reverse strand (2)
   #col4: RC of inner primer sequence

cat ./info/2023-11-30.locus_info2 | awk '{print $3}' | sort | uniq
    #all on the reverse strand this time.





#####7-2, Probe inforamtion (updated)
 # ProbeはPrimerと同じ方向性で長めに設計
 # SNPの片側10bp、もう片側を30bpの計31bpでdesignしていて
 # primer designのinputファイルの3, 4カラム目が2 30の場合には"10-SNP-30"
 # 1 30の場合には"30-SNP-10"という形でdesignしています。

 #THIS JOB needs the output of 8-2 (splitted bam files)

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

ODIR=analysis/2023-11-30/probe
mkdir -p $ODIR

locus_info=info/2023-11-30.locus_info2
sample_barcode_file=info/2023-11-30.sample_barcode_pair

cat $locus_info | awk '{print $1}' | while read snpid;do
   echo $snpid ...
   
   OFILE=$ODIR/$snpid.probe_design.info
   echo "ID Frq Up_seq Allele Dw_seq" > $OFILE
   
   #adjust as needed
   chr=$( echo $snpid | cut -d "_" -f 1 ) #include "chr"
   snppos=$( echo $snpid | cut -d "_" -f 2 )
   refallele=$( echo $snpid | cut -d "_" -f 3 )
   altallele=$( echo $snpid | cut -d "_" -f 4 )
   strand=$( cat $locus_info |
      awk -v snpid=$snpid '{if($1==snpid){print $3}}' ) #the strand of the reverse primer
   win_fow=15 #the opposite side of the reverse primer
   win_rev=5 #the same side of the reverse primer (in this project need to be very small)
   coordinate=GRCh38 #GRCh37
   FASTA=/home/imgishi/reference/Fasta/$coordinate/$coordinate.primary_assembly.genome.bg.fa.gz
   
   #probe start and end position
   if [ $strand = 1 ] ; then 
      pstart=$( expr $snppos - $win_rev )
      pend=$( expr $snppos + $win_fow )
   else
      pstart=$( expr $snppos - $win_fow )
      pend=$( expr $snppos + $win_rev )
   fi
   
   ###1, standard probe using reference sequence
   #Upstream sequence
   from=$pstart
   to=$( expr $snppos - 1 )
   region=$( echo $chr":"$from"-"$to )
   seq_up=$( samtools faidx $FASTA $region |
      perl -pe "s/\n/\t/g" | awk '{print $2}' )
   
   #Downstream sequence
   from=$( expr $snppos + 1 )
   to=$pend
   region=$( echo $chr":"$from"-"$to )
   seq_dw=$(  samtools faidx $FASTA $region |
      perl -pe "s/\n/\t/g" | awk '{print $2}' )
   
   #out
   echo "REF" NA $seq_up $refallele $seq_dw >> $OFILE
   echo "ALT" NA $seq_up $altallele $seq_dw >> $OFILE
   
   ###2, get seq info from bam
   #R2 and read mapped in proper pair (130)
   for jobid in $(seq 1 3);do
      id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
      barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
      BAM=bowtie/20231130/$id.$barcode/$snpid/$snpid.100bp_q30.bam
      
      samtools view -f 130 $BAM |
      awk -v snppos=$snppos -v pstart=$pstart -v pend=$pend '{
         r_start=$4;
         sequence=$10;
         r_length=length(sequence);
         r_end=r_start + r_length - 1;
         if( r_start < pstart && pend < r_end ){
            #up sequence
            r_start_mod = pstart - r_start + 1;
            extract_len = snppos - pstart;
            seq_up = substr(sequence, r_start_mod, extract_len);
            #allele
            r_start_mod = snppos - r_start + 1;
            extract_len = 1;
            allele = substr(sequence, r_start_mod, extract_len);
            #dw sequence
            r_start_mod = snppos - r_start + 2;
            extract_len = pend - snppos;
            seq_dw = substr(sequence, r_start_mod, extract_len);
            print seq_up "_" allele "_" seq_dw
         }
      }' |
      sort | uniq -c | sort -rn |
      head -n 10 |
      awk -v id=$id '{
        split($2,D,"_");
        print id, $1,D[1],D[2],D[3]
      }' >> $OFILE
      
   done
   
done

 #summary file using jupyter
 #analysis/2023-11-30/probe/snp_probe_fw.info
 #SNP ID ALLELE PROBE















#####8-1, Bam split into each region (updated)

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/split_bam-2023-11-30.log \
   -e log/split_bam-2023-11-30.error \
   ./script/split_bam-2023-11-30.sh
    #Your job-array 84418901.1-3:1 ("split_bam-2023-11-30.sh") has been submitted

 #output check:
locus_info=info/2023-11-30.locus_info2
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   cat $locus_info | awk '{print $1}' | while read snpid;do
      if [ ! -s bowtie/20231130/$id.$barcode/$snpid/$snpid.100bp.bam ] ; then echo $id $barcode $snpid ; fi
   done
done

 # ./script/split_bam-2023-11-30.sh
 #-------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
sample_barcode_file=info/2023-11-30.sample_barcode_pair
locus_info=info/2023-11-30.locus_info2
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
ODIR=bowtie/20231130/$id.$barcode
BAM=bowtie/20231130/$id.$barcode/d1.sorted.bam
 #####################

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

cat $locus_info | awk '{print $1}' | while read snpid;do
   CHR=$( echo $snpid | cut -d "_" -f 1 )
   POS=$( echo $snpid | cut -d "_" -f 2 )
   FROM=$( expr $POS - 100 )
   TO=$( expr $POS + 100 )
   if [ $FROM -lt 1 ] ; then FROM=1 ; fi
   Region=$CHR:${FROM}-${TO}
   
   mkdir -p $ODIR/$snpid
   
   samtools view -b $BAM $Region > $ODIR/$snpid/$snpid.100bp.bam
   
   samtools index $ODIR/$snpid/$snpid.100bp.bam
   
done

 #-------------------------------------#


 #memo
snpid=chr1_16613699_G_A
Region=chr1:16613599-16613799
samtools view $BAM $Region |
cut -f 4 | sort -n | head -n 1 #16613492 (probably, end base of reads overlap with the target region)

samtools view $BAM $Region |
cut -f 4 | sort -rn | head -n 1 #16613796 (start base of reads overlap with the target region)

Region=chr1:16613499-16613899
samtools view $BAM $Region |
cut -f 4 | sort -n | head -n 1 #16613393 (probably, end base of reads overlap with the target region)

samtools view $BAM $Region |
cut -f 4 | sort -rn | head -n 1 #16613889 (start base of reads overlap with the target region)









#####8-2, Bam split into each region (updated)
 #Q30 bam file

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/split_bam_q30-2023-11-30.log \
   -e log/split_bam_q30-2023-11-30.error \
   ./script/split_bam_q30-2023-11-30.sh
    #Your job-array 84419745.1-3:1 ("split_bam_q30-2023-11-30.sh") has been submitted

 #output check:
locus_info=info/2023-11-30.locus_info2
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   cat $locus_info | awk '{print $1}' | while read snpid;do
      if [ ! -s bowtie/20231130/$id.$barcode/$snpid/$snpid.100bp_q30.bam ] ; then echo $id $barcode $snpid ; fi
   done
done

 # ./script/split_bam_q30-2023-11-30.sh
 #-------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
sample_barcode_file=info/2023-11-30.sample_barcode_pair
locus_info=info/2023-11-30.locus_info2
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
ODIR=bowtie/20231130/$id.$barcode
BAM=bowtie/20231130/$id.$barcode/d1.sorted.q30.bam
 #####################

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

cat $locus_info | awk '{print $1}' | while read snpid;do
   CHR=$( echo $snpid | cut -d "_" -f 1 )
   POS=$( echo $snpid | cut -d "_" -f 2 )
   FROM=$( expr $POS - 100 )
   TO=$( expr $POS + 100 )
   if [ $FROM -lt 1 ] ; then FROM=1 ; fi
   Region=$CHR:${FROM}-${TO}
   
   mkdir -p $ODIR/$snpid
   
   samtools view -b $BAM $Region > $ODIR/$snpid/$snpid.100bp_q30.bam
   
   samtools index $ODIR/$snpid/$snpid.100bp_q30.bam
   
done

 #-------------------------------------#







#####9-1, QC (updated): Coverage

 ###1, read count
 #R2 and read mapped in proper pair (130)

OFILE=analysis/2023-11-30/read_count/split_bam_r2_count.txt
echo "id barcode snpid nread" > $OFILE

for jobid in $(seq 1 3);do
   #### MODIFY HERE ###
   sample_barcode_file=info/2023-11-30.sample_barcode_pair
   locus_info=info/2023-11-30.locus_info2
   id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
   barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
   ODIR=bowtie/20231130/$id.$barcode
   
   cat $locus_info | awk '{print $1}' | while read snpid;do
      BAM=$ODIR/$snpid/$snpid.100bp_q30.bam
      nread=$( samtools view -f 130 $BAM | wc -l )
      echo $id $barcode $snpid $nread >> $OFILE
   done
done


 ###2, genome wide coverage
 #make bin information
R
options(scipen=999) 

d1 <- read.table("/home/imgishi/reference/STAR_ref/Grch38_G26/chrNameLength.txt")
OFILE="info/grch38_1e7_bin.txt"
bin_size=1e7
bed <- data.frame()
for(chr in 1:22){
   clen <- subset(d1,V1==paste0("chr",chr))$V2
   nbin <- clen / bin_size + 1
   for( i in 1:nbin){
      from=0 + bin_size * ( i - 1 )
      to=from + bin_size
      if(to >= clen){ to=clen}
      from = as.character(from)
      to = as.character(to)
      dump <- data.frame(chr=paste0("chr",chr),from,to)
      bed <- rbind(bed,dump)
   }
}
bed$ID <- paste0("D:",1:nrow(bed))
write.table(bed,OFILE,col.names=F,row.names=F,append=F,quote=F,sep="\t")
    #R off

 #calculate coverage
for jobid in $(seq 1 3);do
   #### MODIFY HERE ###
   sample_barcode_file=info/2023-11-30.sample_barcode_pair
   locus_info=info/2023-11-30.locus_info2
   id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
   barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
   ODIR=bowtie/20231130/$id.$barcode
   BAM=bowtie/20231130/$id.$barcode/d1.sorted.q30.bam
   
   samtools view -b -f 130 $BAM chr{1..22} |
   bedtools coverage -a info/grch38_1e7_bin.txt -b - > analysis/2023-11-30/$id.$barcode.cov
   
done

 #snp and bin
cat info/2023-11-30.locus_info2 |
awk 'BEGIN{OFS="\t"}{
   split($1,D,"_");
   print D[1], D[2] - 1, D[2], $1
}' |
sort -k 1,1 -k 2,2n |
bedtools coverage -a info/grch38_1e7_bin.txt -b - > analysis/2023-11-30/target_snp.cov









#####9-2, QC (updated): insert size
 #will be analyzed in the final summary file
 #$ODIR/$snpid/$snpid.ref_umi_info.*.txt
 #$ODIR/$snpid/$snpid.alt_umi_info.*.txt





#####9-3, QC (updated): Indel check
 #MAPQ > 10
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

sample_barcode_file=info/2023-11-30.sample_barcode_pair
locus_info=info/2023-11-30.locus_info2
ODIR=analysis/2023-11-30/cigar
mkdir -p $ODIR

cat $locus_info | awk '{print $1}' | while read snpid;do
   echo $snpid...
   
   OFILE=analysis/2023-11-30/cigar/$snpid.character_cigar
   
   rm -f $OFILE.tmp
   
   cat $sample_barcode_file | while read line;do
      id=$( echo $line | awk '{print $1}' )
      barcode=$( echo $line | awk '{print $2}' )
      BAM=bowtie/20231130/$id.$barcode/$snpid/$snpid.100bp.bam
      
      #R2, uniquly mapped read, cigar
      #extract only character
      samtools view -F 64 $BAM |
      grep -F -v "XS:i:" - |
      awk '{print $6}' |
      awk '{gsub(/[0-9]/, "",$1); print $1}'  >> $OFILE.tmp
      
   done
   
   cat $OFILE.tmp | sort | uniq -c | sort -rn | awk '{print $2,$1}' > $OFILE
   
   rm -f $OFILE.tmp
   
done

 #overview
R
snplist <- as.character(read.table("info/2023-11-30.locus_info2")[,1])
res <- data.frame()
for( i in 1:length(snplist) ){
   d1 <- read.table(paste0("analysis/2023-11-30/cigar/",snplist[ i ],".character_cigar"))
   colnames(d1) <- c("Cigar","Readcount")
   M_count <- subset(d1,Cigar=="M")$Readcount
   d2 <- subset(d1, Cigar != "M" & Readcount > M_count*0.05 )
   if(nrow(d2)>0){
      d2$M_count <- M_count
      d2$Ratio <- d2$Readcount / d2$M_count
      d2$SNP=snplist[ i ]
      d2 <- d2[,c("SNP","Cigar","Ratio","Readcount","M_count")]
      res <- rbind(res,d2)
   }
}

show(res)
    #             SNP                  Cigar      Ratio          Readcount   M_count
    # chr1_22025454_G_T     MIM     0.07351518     14357    195293
    # chr17_46194146_G_C   MIM    0.10408269     32405     311339
    # chr2_231514523_G_T   MIM     0.49520302    121298    244946
    # chr22_41544436_A_G   MIM     0.40609382     27789     68430
    # chr7_139341719_G_T   MDM    0.74344736     80867    108773
    #these loci might produce incorrect allelic imbalance estimate






#####9-4, QC (updated): Indel check
 #MAPQ > 30
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

sample_barcode_file=info/2023-11-30.sample_barcode_pair
locus_info=info/2023-11-30.locus_info2
ODIR=analysis/2023-11-30/cigar
mkdir -p $ODIR

cat $locus_info | awk '{print $1}' | while read snpid;do
   echo $snpid...
   
   OFILE=analysis/2023-11-30/cigar/$snpid.character_cigar_q30
   
   rm -f $OFILE.tmp
   
   cat $sample_barcode_file | while read line;do
      id=$( echo $line | awk '{print $1}' )
      barcode=$( echo $line | awk '{print $2}' )
      BAM=bowtie/20231130/$id.$barcode/$snpid/$snpid.100bp_q30.bam
      
      #R2, uniquly mapped read, cigar
      #extract only character
      samtools view -F 64 $BAM |
      grep -F -v "XS:i:" - |
      awk '{print $6}' |
      awk '{gsub(/[0-9]/, "",$1); print $1}'  >> $OFILE.tmp
      
   done
   
   cat $OFILE.tmp | sort | uniq -c | sort -rn | awk '{print $2,$1}' > $OFILE
   
   rm -f $OFILE.tmp
   
done

 #overview
R
snplist <- as.character(read.table("info/2023-11-30.locus_info2")[,1])
res <- data.frame()
for( i in 1:length(snplist) ){
   d1 <- read.table(paste0("analysis/2023-11-30/cigar/",snplist[ i ],".character_cigar_q30"))
   colnames(d1) <- c("Cigar","Readcount")
   M_count <- subset(d1,Cigar=="M")$Readcount
   d2 <- subset(d1, Cigar != "M" & Readcount > M_count*0.05 )
   if(nrow(d2)>0){
      d2$M_count <- M_count
      d2$Ratio <- d2$Readcount / d2$M_count
      d2$SNP=snplist[ i ]
      d2 <- d2[,c("SNP","Cigar","Ratio","Readcount","M_count")]
      res <- rbind(res,d2)
   }
}

show(res)
    #             SNP                  Cigar      Ratio          Readcount   M_count
    # chr2_231514523_G_T   MIM     0.4816956    111882       232267
    # chr7_139341719_G_T   MDM   0.7359172     79652        108235
    #much improved compared with Q > 10
    #use Q > 30






#####9-5, QC (updated): Indel check
 #using Hatano bam
 #can we detect loci with many indel before conducting UNIChro-seq?

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH

sample_barcode_file=info/2023-11-30.sample_barcode_pair
locus_info=info/2023-11-30.locus_info2

ODIR=analysis/2023-11-30/cigar_hatano_bam
mkdir -p $ODIR

BAM=/home/ha7477/share/to_imgkono2/result/Jurkat/eightsamplesmerge/eightmerged_markdup_F1804f2q30_autosomal_bfilt_addRG.bam

cat $locus_info | awk '{print $1}' | while read snpid;do
   CHR=$( echo $snpid | cut -d "_" -f 1 )
   POS=$( echo $snpid | cut -d "_" -f 2 )
   FROM=$( expr $POS - 100 )
   TO=$( expr $POS + 100 )
   if [ $FROM -lt 1 ] ; then FROM=1 ; fi
   Region=$CHR:${FROM}-${TO}
   
   TMPBAM=tmp.100bp.bam
   samtools view -b $BAM $Region > $TMPBAM
   
   OFILE=$ODIR/$snpid.character_cigar
   
   samtools view -F 64 $TMPBAM |
   grep -F -v "XS:i:" - |
   awk '{print $6}' |
   awk '{gsub(/[0-9]/, "",$1); print $1}' |
   sort | uniq -c | sort -rn | awk '{print $2,$1}' > $OFILE
   
   rm -f $TMPBAM
   
done

 #overview
R
snplist <- as.character(read.table("info/2023-11-30.locus_info2")[,1])
res <- data.frame()
for( i in 1:length(snplist) ){
   d1 <- read.table(paste0("analysis/2023-11-30/cigar_hatano_bam/",snplist[ i ],".character_cigar"))
   colnames(d1) <- c("Cigar","Readcount")
   M_count <- subset(d1,Cigar=="M")$Readcount
   d2 <- subset(d1, Cigar != "M" & Readcount > M_count*0.05 )
   if(nrow(d2)>0){
      d2$M_count <- M_count
      d2$Ratio <- d2$Readcount / d2$M_count
      d2$SNP=snplist[ i ]
      d2 <- d2[,c("SNP","Cigar","Ratio","Readcount","M_count")]
      res <- rbind(res,d2)
   }
}

show(res)
     #              SNP                  Cigar      Ratio         Readcount  M_count
     # chr1_16613699_G_A    MS         0.07692308         3       39
     #  chr1_22025454_G_T    MS        0.14912281        17      114
     # chr1_22025454_G_T     SM        0.07017544         8       114
     # chr17_46194146_G_C   MIM      0.08095238        17       210
     # chr17_46194146_G_C    SM       0.05714286        12       210
     # chr17_57551724_G_A    MS        0.05645161         7       124
     # chr18_45687443_T_G    MDM     0.28735632        25       87
     # chr2_231514523_G_T   MIM        0.15662651        26      166
     # chr2_231514523_G_T    SM         0.05421687         9       166
     # chr22_41544436_A_G    SM         0.05617978         5       89
     # chr22_41544436_A_G    MIM       0.05617978         5        89
     # chr4_118278629_G_C    SM         0.06321839        11      174
     # chr7_139341719_G_T    MDM      0.22222222        28      126
     #substantially different...
     #probably, due to bowtie2 parameter differences




'GCCGCGCAGGCGGAG'
'GCCGCGCAGGCAGAG'





#####10-1, Extract data required for UMI counting (updated)

qsub -pe def_slot 1 \
   -l s_vmem=10G,mem_req=10G \
   -cwd \
   -t 1:3 -tc 100 \
   -o log/umi_count_prep-2023-11-30.log \
   -e log/umi_count_prep-2023-11-30.error \
   ./script/umi_count_prep-2023-11-30.sh
    #Your job-array 84429439.1-3:1 ("umi_count_prep-2023-11-30.sh") has been submitted

 #output check: all 11, ok!
locus_info=info/2023-11-30.locus_info2
cat info/2023-11-30.sample_barcode_pair | while read line;do
   id=$( echo $line | awk '{print $1}' )
   barcode=$( echo $line | awk '{print $2}' )
   cat $locus_info | awk '{print $1}' | while read snpid;do
      ll bowtie/20231130/$id.$barcode/$snpid/$snpid.ref_umi_info.ED*.txt.gz | wc -l
      ll bowtie/20231130/$id.$barcode/$snpid/$snpid.alt_umi_info.ED*.txt.gz | wc -l
   done
done

 # ./script/umi_count_prep-2023-11-30.sh
 #-------------------------------------#
#!/bin/sh
jobid=${SGE_TASK_ID}

 #### MODIFY HERE ###
sample_barcode_file=info/2023-11-30.sample_barcode_pair
locus_info=info/2023-11-30.locus_info2
id=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $2}}' )
barcode=$( cat -n $sample_barcode_file | awk -v jobid=$jobid '{if($1==jobid){print $3}}' )
ODIR=bowtie/20231130/$id.$barcode
 #####################

ORGPATH=$( echo $PATH )
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

cat $locus_info | awk '{print $1}' | while read snpid;do
   BAM=$ODIR/$snpid/$snpid.100bp_q30.bam
   
   ref_probe=$( cat analysis/2023-11-30/probe/snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="REF" ){print $3}}' )
   alt_probe=$( cat analysis/2023-11-30/probe/snp_probe_fw.info |
      awk -v snpid=$snpid '{if( $1==snpid && $2=="ALT" ){print $3}}' )
   
   #Strategy:
   #Focus on the R1 and R2 1st base map position
   # - R2 should be same as primer site
   # - when R2 is on the revserse strand, R2 1st base = end position in the bam file
   # - when R2 is on the foward strand, R2 1st base = start position in the bam file
   #Flag 146: R2, on REVERSE strand, read mapped in proper pair
   #Flag 162: R2, on FOWARD strand, read mapped in proper pair
   #exclude multiple mapping (R2: -F 64; uniquely mapped: not "XS:i:")
   #prepare input file with basic map data
   
   export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH
   
   ### REF allele ###
   samtools view -f 146 $BAM |
   grep -F -v "XS:i:" - |
   awk 'BEGIN{FS="\t"}{
      Strand="-";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R1_1st_base=$8;
      size=$9;
      if(size< 0){size = -size};
      R2_1st_base = R1_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
   }' |
   grep -F $ref_probe > $ODIR/$snpid/$snpid.ref_umi_info.tmp
   
   samtools view -f 162 $BAM |
   grep -F -v "XS:i:" - |
   awk 'BEGIN{FS="\t"}{
      Strand="+";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R2_1st_base=$4;
      size=$9;
      if(size< 0){size = -size};
      R1_1st_base = R2_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
   }' |
   grep -F $ref_probe >> $ODIR/$snpid/$snpid.ref_umi_info.tmp
   
   ### ALT allele ###
   samtools view -f 146 $BAM |
   grep -F -v "XS:i:" - |
   awk 'BEGIN{FS="\t"}{
      Strand="-";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R1_1st_base=$8;
      size=$9;
      if(size< 0){size = -size};
      R2_1st_base = R1_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
   }' |
   grep -F $alt_probe > $ODIR/$snpid/$snpid.alt_umi_info.tmp
   
   samtools view -f 162 $BAM |
   grep -F -v "XS:i:" - |
   awk 'BEGIN{FS="\t"}{
      Strand="+";
      readid=$1;
      split(readid,D,"_");
      UMI=D[2];
      CIGAR=$6;
      gsub(/[0-9]/, "",CIGAR);
      R2_1st_base=$4;
      size=$9;
      if(size< 0){size = -size};
      R1_1st_base = R2_1st_base + size - 1;
      sequence=$10;
      print UMI, R1_1st_base, R2_1st_base, CIGAR, size, sequence, Strand
   }' |
   grep -F $alt_probe >> $ODIR/$snpid/$snpid.alt_umi_info.tmp
   
   #R script
   #remove redundant information using UMI
   #input: $ODIR/$snpid/$snpid.ref_umi_info.tmp
   #outout: $ODIR/$snpid/$snpid.ref_umi_info.ED0-10
   export PATH=/home/imgishi/miniconda3/envs/r4.2.2/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH
   
   Nrow=$( cat $ODIR/$snpid/$snpid.ref_umi_info.tmp | wc -l )
   if [ $Nrow -gt 1 ] ; then
      R --vanilla --slave --args $ODIR/$snpid/$snpid.ref_umi_info.tmp \
      < /home/imgishi/wd/crisprqtl/script/umi_count_v1.R
   else
      echo "No input in $ODIR/$snpid/$snpid.ref_umi_info.tmp"
   fi

   Nrow=$( cat $ODIR/$snpid/$snpid.alt_umi_info.tmp | wc -l )
   if [ $Nrow -gt 1 ] ; then
      R --vanilla --slave --args $ODIR/$snpid/$snpid.alt_umi_info.tmp \
      < /home/imgishi/wd/crisprqtl/script/umi_count_v1.R
   else
      echo "No input in $ODIR/$snpid/$snpid.alt_umi_info.tmp"
   fi
   
   rm -f $ODIR/$snpid/$snpid.ref_umi_info.tmp
   rm -f $ODIR/$snpid/$snpid.alt_umi_info.tmp
   
   gzip -f $ODIR/$snpid/$snpid.ref_umi_info.*.txt
   gzip -f $ODIR/$snpid/$snpid.alt_umi_info.*.txt
   
done

 #-------------------------------------#

 # /home/imgishi/wd/crisprqtl/script/umi_count_v1.R
 #-------------------------------------#
args <- commandArgs(trailingOnly = T)
IFILE <- as.character(args[1])
umi_length=17 #modify as needed
 # IFILE="bowtie/20231130/01_S1_L001.CTCTCTAT/chr1_16613699_G_A/chr1_16613699_G_A.alt_umi_info.tmp"
OFILE_header <- gsub(".tmp$","",IFILE)

library(data.table)
library(magrittr)

 #read in data
d1 <- fread(IFILE)
d1 <- as.data.frame(d1)
colnames(d1) <- c("UMI","R1_1st_base","R2_1st_base","CIGAR","Size","R2_Seq","R2_strand")
    #nrow(d1) #62424
    #head(d1,n=3)

 #add paired position tag
d1$ptag <- paste0(d1$UMI,":",d1$R1_1st_base,":",d1$R2_1st_base)

 #add tag frequency data
TB <- table(d1$ptag)
Frq <- data.frame(ptag=names(TB),ptag_count=c(TB))
d2 <- merge(d1,Frq,by="ptag")
    #nrow(d2) #62424
if( nrow(d1) != nrow(d2) ){print("ERROR in counting reads")}

 #UMI data frame
TB <- table(d1$UMI)
umi_df <- data.frame(umi=names(TB),count=c(TB))

 #QC1
is_valid_umi <- sapply(umi_df$umi,function(x){
    x <- gsub("A","",x);
    x <- gsub("T","",x);
    x <- gsub("G","",x);
    x <- gsub("C","",x);
    if(nchar(x)==0){return(TRUE)}else{return(FALSE)}
})
    #nrow(umi_df) #15067
umi_df <- umi_df[is_valid_umi, ]
    #nrow(umi_df) #15067

 #QC2
is_valid_umi <- sapply(umi_df$umi,function(x){
    if(nchar(x)==umi_length){return(TRUE)}else{return(FALSE)}
})
    #nrow(umi_df) #15067
umi_df <- umi_df[is_valid_umi, ]
    #nrow(umi_df) #15067

 #STEP0
umi_df <- umi_df[order(umi_df$count, decreasing=TRUE), ]
    #head(umi_df,n=3)

 #STEP1: Convert ATGC into numeric form
umi_df$umi_digit <- sapply(umi_df$umi,function(x){
    x <- gsub("A","1000",x);
    x <- gsub("T","0100",x);
    x <- gsub("G","0010",x);
    x <- gsub("C","0001",x);
    return(x)
})
    #head(umi_df$umi_digit,n=3)

umi_matrix <- sapply(umi_df$umi_digit,function(x){
    as.numeric(unlist(strsplit(split="",x)))
})
colnames(umi_matrix) <- as.character(umi_df$umi)
    #dim(umi_matrix) #row: 4 x nchar(UMI)=68, column: length(UMI)
    #umi_matrix[1:4,1:4]

 #STEP2: Hamming distance matrix 
n_bases_matched <- t(umi_matrix) %*% umi_matrix
n_bases_mismatched <- umi_length - n_bases_matched
    #dim(n_bases_mismatched) #length(UMI) x length(UMI)
    #n_bases_mismatched[1:4,1:4]

 ### UMI counting with no error correction ###
edit_distance=0
res <- data.frame(
   index_umi=umi_df$umi,
   count_sum=umi_df$count
)

 #summary: filter by index umi, unique UMI-R2pos pairs
d3 <- d2[ d2$UMI %in% res$index_umi, ]
d4 <- d3[!duplicated(d3$ptag),]
d5 <- merge(res,d4,by.x="index_umi",by.y="UMI")
      # nrow(d4)==nrow(d5) #TRUE!
OFILE=paste0(OFILE_header,".ED",edit_distance,".txt")
write.table(d5, OFILE, col.names=T, row.names=F, append=F, quote=F, sep="\t")

 ### UMI counting with error correction ###
for( edit_distance in 1:10 ){
   umi_df_mod <- umi_df
   res <- data.frame()
   
   for(i in 1:nrow(umi_df)){
       if( nrow(umi_df_mod)==0 ){ 
           break 
       }else{
           #index umi
           index_umi <- umi_df_mod$umi[1]
           
           #hamming distance from index umi
           diff_with_index <- n_bases_mismatched[index_umi,]
           
           #all umis probably coming from the index umi (index umi cluster)
           index_umi_cluster <- names(diff_with_index[diff_with_index <= edit_distance])
           
           #sum of count of the index umi cluster
           count_sum <- umi_df[umi_df$umi %in% index_umi_cluster, "count"] %>% sum()
           
           #umi_df_mod for the next cycle (remove index_umi_cluster)
           umi_df_mod <- umi_df_mod[ ! umi_df_mod$umi %in% index_umi_cluster, ]
           
           #output
           dump <- data.frame(index_umi,count_sum)
           res <- rbind(res,dump)
       }
   }
   
   #summary: filter by index umi, unique UMI-R2pos pairs
   d3 <- d2[ d2$UMI %in% res$index_umi, ]
   d4 <- d3[!duplicated(d3$ptag),]
   d5 <- merge(res,d4,by.x="index_umi",by.y="UMI")
      # nrow(d4)==nrow(d5) #TRUE!
   OFILE=paste0(OFILE_header,".ED",edit_distance,".txt")
   write.table(d5, OFILE, col.names=T, row.names=F, append=F, quote=F, sep="\t")
}

 #-------------------------------------#







#####10-2, UMI counting (updated)
 #using jupyter 
 #for ulimit in 100 150 200 250 300 400 500;do



