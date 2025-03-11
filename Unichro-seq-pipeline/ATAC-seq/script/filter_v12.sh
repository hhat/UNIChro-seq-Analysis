#!/bin/sh

# Parameters for filtering
size_cutoff=${1}
count_cutoff=${2}
batch=${3}

# Environment setting
ORGPATH=$( echo $PATH )
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1
export PATH=/home/imgishi/miniconda3/envs/genetics1/bin:/home/imgishi/miniconda3/condabin:/home/imgishi/miniconda3/bin:/home/imgishi/wd/crisprqtl/script:$ORGPATH
   
# Working directory
wd=/home/imgkono/wd/img/crispr_qtl/UNIChro_seq2/

# Find all UMI data files for the specified batch
files=`ls ${wd}/bowtie/${batch}/*/chr*/*_umi_info.ED2_index_umi_final.txt`

# Process each file
for input_file in $files
do
  # Define output filename with filter parameters
  output_file=${input_file}_size${size_cutoff}_count${count_cutoff}_halffilter.txt
  
  # Skip if output file already exists
  if [ -f "$output_file" ]; then
    echo "Output file $output_file already exists. Skip"
    continue
  fi
  
  # Check if input file exists
  if [ ! -f "$input_file" ]; then
    echo "File $input_file does not exist. Skip"
    continue
  fi
  
  # Apply filtering criteria using awk
  cat "$input_file" | \
  awk -v size_cutoff="$size_cutoff" -v count_cutoff="$count_cutoff" '
  BEGIN {
    FS = "\t" # Set field separator to tab
  }
  NR==1 {
    # Get column positions from header row
    for (i=1; i<=NF; i++) {
      if ($i == "Size_mode") size_col = i
      if ($i == "count") count_col = i
      if ($i == "refaltchanged_count") refaltchanged_count_col = i
      if ($i == "samefastq_bg_woself_count") samefastq_bg_woself_count_col = i
      if ($i == "umi_count_barcode") umi_count_barcode_col = i
    }
    print # Output header row
  }
  NR>1 {
    # Apply filtering criteria using column positions
    if ($size_col <= size_cutoff && 
        $count_col >= count_cutoff && 
        $refaltchanged_count_col < $count_col && 
        $samefastq_bg_woself_count_col < $count_col && 
        $umi_count_barcode_col < $count_col) {
      print
    }
  }' > "$output_file"
done