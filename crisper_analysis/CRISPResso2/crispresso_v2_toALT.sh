#!/bin/sh
#
# CRISPResso2 analysis script for allele-specific editing
# This is an example script for toALT. For toREF analysis, follow the same approach.

snp_line=${SGE_TASK_ID}

# Environment setup
export PATH=/home/ha7477/tools/miniconda3/bin:${PATH}
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1
wd=/home/ha7477/works/umi/crispresso
gRNA_file=/home/ha7477/works/umi/crispresso/data/R10_gRNA.txt
prob_file=/home/imgkono/wd/img/crispr_qtl/DNAseq2/info/rR10_snp_probe_fw.info 
eval "$(/home/ha7477/tools/miniconda3/bin/conda shell.bash hook)"
conda activate crispresso2


{
  tail -n +2 /home/ha7477/works/umi/info/v1_DNAseq.txt
  tail -n +2 /home/ha7477/works/umi/info/v2_DNAseq.txt
} | while IFS=$'\t' read -r sample_raw replicate sample_id barcode bam_day; do
    # Format sample identifier
    rep_num=$(printf "%02d" $replicate)
    sample_id2="${sample_raw}_REP${rep_num}"
    
    echo "Processing sample_id: $sample_id, sample_id2: $sample_id2, bam_day: $bam_day"

    # Skip if not type C or R
    type=$(echo $sample_raw | cut -d'_' -f3)
    if [ "$type" = "A" ] || [ "$type" != "C" ] && [ "$type" != "R" ]; then
        continue
    fi
    echo "Type is $type. Processing..."
    
    # Process the current SNP
    tail -n +2 "${gRNA_file}" | head -n ${snp_line} | tail -n 1 | while IFS=$'\t' read -r rezaid hg19pos hg38pos rsid jurkat guide strand; do
        # Reset variables
        AMPLICON_SEQUENCE=""
        EXPECTED_SEQUENCE=""

        # Convert hg38 position to SNP format
        snp=$(echo "$hg38pos" | sed 's/^/chr/' | tr ':' '_')
        snp2="$rsid"
        GUIDE_SEQUENCE="$guide"
        od=${wd}/result_afterprobeQC/${snp2}/${sample_id2}
        mkdir -p $od
        cd ${od}
       
        # Reverse complement if on negative strand
        GUIDE_SEQUENCE2=`echo $GUIDE_SEQUENCE | awk 'BEGIN {
            comp["A"] = "T"; comp["T"] = "A"; comp["G"] = "C"; comp["C"] = "G";
        }
        {
            seq = $0; rev = "";
            for (i = length(seq); i > 0; i--) {
                base = substr(seq, i, 1);
                rev = rev comp[base];
            }
            print rev;
        }'`

        if [ "$strand" = "-" ]; then
            GUIDE_SEQUENCE="$GUIDE_SEQUENCE2"
        fi   

        echo "Processing ${snp} (${snp2})..."
        
        # Get BAM file path
        bam="/home/imgkono/wd/img/crispr_qtl/DNAseq2/bowtie/${bam_day}/${sample_id}/${snp}/${snp}.R2_after_probe.bam"
     
        # Extract most frequent sequences
        samtools view ${bam} | cut -f 10 | sort | uniq -c | sort -nr | head -n 5 > temp_sequences.txt
        
        # Get reference and alternate sequences
        ref_seq=$(awk -v snp="$snp" '$1==snp && $2=="REF" {print $3}' "$prob_file")
        alt_seq=$(awk -v snp="$snp" '$1==snp && $2=="ALT" {print $3}' "$prob_file")

        # Find amplicon sequence containing reference
        while read count seq; do
            if echo "$seq" | grep -q "$ref_seq"; then
                AMPLICON_SEQUENCE="$seq"
                break
            fi
        done < temp_sequences.txt
        
        # Generate expected sequence with ALT
        EXPECTED_SEQUENCE=${AMPLICON_SEQUENCE//$ref_seq/$alt_seq}
        
        # Verify guide sequence is in amplicon
        if ! echo "$AMPLICON_SEQUENCE" | grep -q "$GUIDE_SEQUENCE"; then
            echo "Error: GUIDE_SEQUENCE not found in AMPLICON_SEQUENCE for ${snp}"
            continue
        fi
        
        # Calculate quantification window center
        snp_pos=$(echo $snp | cut -d'_' -f2)
        bam_pos=$(samtools view ${bam} | awk -v seq="$AMPLICON_SEQUENCE" '$10==seq {print $4}' | sort | uniq -c | sort -nr | head -n 1 | awk '{print $2}')
        pos_diff=$((snp_pos - bam_pos--))
        guide_pos=$(echo "$AMPLICON_SEQUENCE" | grep -b -o "$GUIDE_SEQUENCE" | cut -d: -f1)
        guide_length=${#GUIDE_SEQUENCE}
        CENTER=$((pos_diff - (guide_pos + guide_length -1)))
        
        echo "SNP position: $snp_pos, BAM start: $bam_pos, Diff: $pos_diff"
        echo "Guide position: $guide_pos, Guide length: $guide_length"
        echo "CENTER: $CENTER"    
            
        # Run CRISPResso
        CRISPResso --bam_input ${bam} \
            --amplicon_seq "${AMPLICON_SEQUENCE}" \
            --expected_hdr_amplicon_seq "${EXPECTED_SEQUENCE}" \
            --guide_seq "${GUIDE_SEQUENCE}" \
            --quantification_window_size 20 \
            --quantification_window_center ${CENTER} \
            -n ${snp2}_toALT \
            -o ${od} 

        rm temp_sequences.txt
    done
done