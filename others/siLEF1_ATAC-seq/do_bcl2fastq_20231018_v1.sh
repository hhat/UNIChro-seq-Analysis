#!/bin/bash
#$ -S /bin/bash

/usr/local/package/bcl2fastq/2.20.0.422/bin/bcl2fastq \
--runfolder-dir /home/imgkawa/data/img/novaseq/20231016_FC01127/231013_A00750_0582_AHTGH3DRX2/ \
--sample-sheet /home/ha7477/share/to_imgkono2/data/FC01127_sample_sheet.csv \
--output-dir /home/imgkawa/data/img/novaseq/20231016_FC01127/fastq/231013_A00750_0582_AHTGH3DRX2/ \
--use-bases-mask y100n,I8,I8,y100n \
--processing-threads 28 \
--loading-threads 4 \
--writing-threads 4
