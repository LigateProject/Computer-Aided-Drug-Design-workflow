#!/bin/bash
set -e

OUTPUT_PATH=${ligate}/CADDValidation
PATH_TO_SCRIPTS=CADD/scripts

targets=(bace_p2)

cd ${OUTPUT_PATH}
for target in ${targets[@]}; do
cd ${target}
for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}

sbatch ${PATH_TO_SCRIPTS}/equilibrate_complex.sh
sbatch ${PATH_TO_SCRIPTS}/equilibrate_ligand.sh

cd ..
done # poses
cd ..
done # edges
cd ..
done # targets
