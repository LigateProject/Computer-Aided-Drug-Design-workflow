#!/bin/bash

#SBATCH --job-name=CADDInputFileGenerationPart2
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH -t 24:00:00
#SBATCH -A cin_staff
#SBATCH -o slurm-%J.out

set -e

# path to additional modules
module use ${CADD_SOFTWARE_MODULES}

OUTPUT_PATH=${CADD_OUTPUT_DIR}
PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

target=${CADD_TARGET}

cd ${OUTPUT_PATH}/${target}

# check result of energy minimisation
for edge in edge_*; do
cd ${edge}
# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the energy minimisation could not be completed successfully for any of the poses. Deleting directory!"
cd ..
rm -rf ${edge}
fi
cd ..
done

# remove directory if none of the edges survives
if (( $(ls -lh | grep "edge_" | wc -l) == 0 ));
then
echo "The energy minimisation could not be completed successfully for any of the edges. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# prepare TPR files for equilibration
bash ${PATH_TO_SCRIPTS}/prepareEquilibration.sh

echo "CADD: Exiting runInputFileGeneration2.sh" >> $LOGFILE
