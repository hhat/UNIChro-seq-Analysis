#!/bin/bash

## Usage: ./workflow.sh <input_file> <output_directory> <parameters>

INPUT_FILE=$1
OUTPUT_DIR=$2
PARAMETERS=${3:-"Tm=62;Window=50;CVG=0.8"}

# script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/script"

mkdir -p $OUTPUT_DIR

cp -f "$INPUT_FILE" "$OUTPUT_DIR/"

INPUT_BASENAME=$(basename $INPUT_FILE)

# step1 creating FASTA file
$SCRIPT_DIR/prepare_fastq.sh $OUTPUT_DIR/$INPUT_BASENAME $OUTPUT_DIR

# step2 design primer
FULL_PARAMS="$PARAMETERS;ODIR=$OUTPUT_DIR;IFILE=$OUTPUT_DIR/$INPUT_BASENAME"
$SCRIPT_DIR/design_primer.sh "$FULL_PARAMS"

# step3: evaluation
$SCRIPT_DIR/evaluate_primer.sh $OUTPUT_DIR

# step4: summary
Rscript $SCRIPT_DIR/make_summary.R $OUTPUT_DIR/$INPUT_BASENAME $OUTPUT_DIR

echo "- Result summary: $OUTPUT_DIR/target.primers.mapinfo_oneline.txt"