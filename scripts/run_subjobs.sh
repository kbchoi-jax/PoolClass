#!/bin/bash
#PBS -l nodes=1:ppn=4,mem=16gb,walltime=23:59:59

cd $PBS_O_WORKDIR

# Add environment modules here, for example ...
module load hdf5/1.8.14
module load R/3.5.1

# Run R script
CASE_ID=`printf %05d $PBS_ARRAYID`
Rscript ${RFILE} _chunk.${CASE_ID} _scrate_elpd_loo.${CASE_ID} ${CORES} ${SEED}
