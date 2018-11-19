#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=32gb,walltime=23:59:00

source activate tenxt
cd $PBS_O_WORKDIR

tenxt extra-zero-test -o ${OUTFILE} ${INFILE}

source deactivate
