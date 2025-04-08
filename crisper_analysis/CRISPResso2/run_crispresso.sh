#!/bin/bash

qsub -pe def_slot 4 -l s_vmem=4G -cwd -t 1:9 crispresso_v2_toALT.sh
