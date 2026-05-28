#!/bin/bash

# Full pipeline script for UNIChro-seq Nested Multiplex Primer Design
# Usage: ./workflow.sh <input_file> <output_directory> <parameters>

INPUT_FILE=$1
OUTPUT_DIR=$2
PARAMETERS=${3:-"Tm=62;Window=50;CVG=0.8"}

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/script"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Copy input file to output directory if not already there
if [ "$INPUT_FILE" != "$OUTPUT_DIR/$(basename $INPUT_FILE)" ]; then
    cp $INPUT_FILE $OUTPUT_DIR/
fi
INPUT_BASENAME=$(basename $INPUT_FILE)

echo "=========================================="
echo "UNIChro-seq Nested Multiplex Primer Design"
echo "=========================================="
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Parameters: $PARAMETERS"
echo "=========================================="

# Input file preparation is completed before running this workflow.
echo "Input data preparation completed"

# STEP 1: Create FASTA files
echo "[STEP 1] Creating FASTA files..."
$SCRIPT_DIR/01_prepare_fastq.sh $OUTPUT_DIR/$INPUT_BASENAME $OUTPUT_DIR

# STEP 2: Design primers
echo "[STEP 2] Designing primers..."
FULL_PARAMS="$PARAMETERS;ODIR=$OUTPUT_DIR;IFILE=$OUTPUT_DIR/$INPUT_BASENAME"
$SCRIPT_DIR/02_design_primer.sh "$FULL_PARAMS"

# STEP 3: Evaluate primer mapping
echo "[STEP 3] Evaluating primer mapping..."
$SCRIPT_DIR/03_evaluate_primer.sh $OUTPUT_DIR

# STEP 4: Summarize results
echo "[STEP 4] Summarizing results..."
Rscript $SCRIPT_DIR/04_make_summary.R $OUTPUT_DIR/$INPUT_BASENAME $OUTPUT_DIR

echo "=========================================="
echo "Pipeline completed"
echo "Results are available in: $OUTPUT_DIR"
echo "- Result summary: $OUTPUT_DIR/target.primers.mapinfo_oneline.txt"
echo "==========================================="
