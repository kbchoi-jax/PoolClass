#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=32gb,walltime=23:59:00

source activate poolclass
cd $PBS_O_WORKDIR

poolclass extra-zero-test --common-scale ${SCALE} -o ${OUTFILE} ${INFILE}

source deactivate
