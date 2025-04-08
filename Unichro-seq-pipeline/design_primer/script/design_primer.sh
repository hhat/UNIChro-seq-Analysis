#!/bin/sh
parameters=$1

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_DOMAIN_NUM_THREADS=1

script_dir=$(dirname "$0")

Rscript ${script_dir}/design_primer.R $parameters
    #"Tm=62;Window=50;CVG=0.8;ODIR=target/demo": parameters (optimise as needed)
    #Tm: melting temperature of primer +/- 3 from this value
    #Dist_snp is not included in this version (specified in the input file)
    #Window: window of site where we design primers
    #CVG: 0-1. required.cvg in openPrimeR. Defines the desired coverage ratio of the templates. 1 indicates it design primers for all sites by relaxing the constrains. 0 indicates it design primers without relaxing the constrains.
    #ODIR: the output directory
    #output file = $ODIR/target.primers
