#!/bin/bash
set -e

OUTPUT_PATH=${CADD_OUTPUT_DIR}
PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

target=${CADD_TARGET}

cd ${OUTPUT_PATH}/${target}
for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}

sbatch --wait ${PATH_TO_SCRIPTS}/energyMinimisation_complex.sh &
sbatch --wait ${PATH_TO_SCRIPTS}/energyMinimisation_ligand.sh &

cd ..
done # poses
cd ..
done # edges
wait
echo "CADD: All energy minimisations finished " >> $LOGFILE

