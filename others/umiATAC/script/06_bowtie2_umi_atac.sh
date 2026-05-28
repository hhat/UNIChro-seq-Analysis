#!/bin/sh

jobid=${SGE_TASK_ID}

wd=${wd:-/home/ha7477/works/umi/umiATAC}
batch=${batch:-20251226}
samplefile=${samplefile:-${wd}/info/2025-12-26.samples}
threads=${threads:-4}
index=${index:-/home/imgishi/reference/Bowtie2_ref/GRCh38/GRCh38}

sample_id=$( cat -n "$samplefile" | awk -v jobid="$jobid" '{if($1==jobid){print $2}}' )
F1=${wd}/cutadapt/${batch}/${sample_id}/${sample_id}.postqc.umi.R1.fastq.gz
F2=${wd}/cutadapt/${batch}/${sample_id}/${sample_id}.postqc.umi.R2.fastq.gz
ODIR=${wd}/bowtie/${batch}/${sample_id}
OFILE=${ODIR}/d1

export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$PATH
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

mkdir -p "$ODIR"

bowtie2 -x "$index" \
  -1 "$F1" \
  -2 "$F2" \
  --threads "$threads" \
  --no-discordant \
  --maxins 1000 \
  -S "$OFILE".sam

samtools view -bS "$OFILE".sam > "$OFILE".bam
samtools view -f 2 -q 10 "$OFILE".bam -o "$OFILE".qc.bam
samtools sort "$OFILE".qc.bam -o "$OFILE".sorted.bam
samtools index "$OFILE".sorted.bam

rm -f "$OFILE".sam "$OFILE".bam "$OFILE".qc.bam

samtools view -f 2 -q 30 "$OFILE".sorted.bam -o "${ODIR}/d1.sorted.q30.bam"
samtools index "${ODIR}/d1.sorted.q30.bam"
