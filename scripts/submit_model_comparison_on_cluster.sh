#!/bin/bash
#PBS -l nodes=1:ppn=4,mem=32gb,walltime=23:59:00

HNAME=`hostname`
if [[ ${HNAME} == *"cadillac"* ]]; then
    module load Anaconda
fi

cd $PBS_O_WORKDIR
source activate scrate

#Rscript ${RFILE} ${LOOMFILE} ${START} ${END} ${OUTBASE}.${START}-${END}.rds
Rscript ${RFILE} ${INFILE} ${OUTFILE} ${CORES} ${SEED}
