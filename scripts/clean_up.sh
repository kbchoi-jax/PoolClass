#!/bin/bash

echo "Cleaning up the following directory: ${WORKDIR}"

echo "Combining error log files..."
cat ${WORKDIR}/submit_model_comparison_on_cluster.sh.e* > ${WORKDIR}/submit_model_comparison_on_cluster.sh.stderr

echo "Combining output log files..."
cat ${WORKDIR}/submit_model_comparison_on_cluster.sh.o* > ${WORKDIR}/submit_model_comparison_on_cluster.sh.stdout

echo "Removing error log files..."
rm ${WORKDIR}/submit_model_comparison_on_cluster.sh.e*

echo "Removing output log files..."
rm ${WORKDIR}/submit_model_comparison_on_cluster.sh.o*

echo "Removing temporary data files..."
rm ${WORKDIR}/_chunk*

echo "Removing temporary result files..."
rm ${WORKDIR}/_scrate*

echo "Done."
