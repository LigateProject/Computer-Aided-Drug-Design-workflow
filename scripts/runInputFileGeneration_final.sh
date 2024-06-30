#!/bin/bash

#SBATCH --job-name=CADDInputFileGeneration
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -o slurm-%J.out

set -e

echo "CADD: Entering runInputFileGeneration_final.sh" >> $LOGFILE
# path to additional modules
module use ${CADD_SOFTWARE_MODULES}

OUTPUT_PATH=${CADD_OUTPUT_DIR}
PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

target=${CADD_TARGET}

export CADD_SUBMISSION_DIR=$(pwd)

# work in directory local to node not to overwhelm file system
mkdir -p ${CADD_LOCAL_DIR}
cd ${CADD_LOCAL_DIR}
if [ -d ${target} ]
then
echo "Another part of the CADD workflow may be working on the same target in the same directory. Stopping this workflow execution not to corrupt the other instance!"
exit 0
fi
cp -r ${OUTPUT_PATH}/${target} ${target}
cd ${target}

bash ${PATH_TO_SCRIPTS}/GROMACS_LiGen_integration_final.sh

# only deal with targets that were not removed by previous stages of the workflow
if ! [ -d ../${target} ]
then
rm -r ${OUTPUT_PATH}/${target}
exit 0
fi

# copy data to OUTPUT_PATH
cd ..
cp -r ${target} ${OUTPUT_PATH}/${target}_new
rm -rf ${OUTPUT_PATH}/${target}
mv ${OUTPUT_PATH}/${target}_new ${OUTPUT_PATH}/${target}
rm -rf ${target}

echo "CADD: Exiting runInputFileGeneration_final.sh" >> $LOGFILE
