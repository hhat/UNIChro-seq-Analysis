#!/bin/bash
#
# CRISPResso2 workflow for CRISPR-Cas9 editing analysis
# Author: Hiroaki Hatano

cd /home/ha7477/works/umi/crispresso

qsub -pe def_slot 4 \
   -l s_vmem=4G \
   -cwd \
   -t 1:9 \
   crispresso_v2_toALT.sh