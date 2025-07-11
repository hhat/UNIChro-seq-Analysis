#!/bin/bash
#$ -S /bin/bash

args=$@

~/tools/miniconda3/envs/de/bin/Rscript --slave --vanilla $args  --max-ppsize=500000
