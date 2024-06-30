#!/bin/bash

#SBATCH --job-name=prepareProductionSimulations
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH --mem=8GB
#SBATCH -t 06:00:00
#SBATCH -o slurm-%J.out

set -e

echo "CADD: PrepareProductionSimulation starting " >>$LOGFILE
module use $CADD_SOFTWARE_MODULES
module load gromacs/2023.2

PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}
PATH_TO_MDP=${CADD_SCRIPTS_DIR}/mdp

target=${CADD_TARGET}

OUTPUT_PATH=${CADD_OUTPUT_DIR}

cd ${OUTPUT_PATH}/${target}

# check equilibration results and clean up first
for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}

## error catching
## The final .gro file of the equilibration must exist
if ! [ -f equiNVT_complex.gro ]
then
echo "For ${pose} of ${edge}, the complex could not be equilibrated. Deleting folder!"
cd ..
rm -rf ${OUTPUT_PATH}/${target}/${edge}/${pose}
continue
fi
if ! [ -f equiNVT_ligand.gro ]
then
echo "For ${pose} of ${edge}, the ligand could not be equilibrated. Deleting folder!"
cd ..
rm -rf ${OUTPUT_PATH}/${target}/${edge}/${pose}
continue
fi

mv equiNVT_complex.gro start_complex.gro
mv equiNVT_ligand.gro start_ligand.gro
rm ions_complex.gro ions_ligand.gro equi* slurm-*

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh ${OUTPUT_PATH}/${target}/${edge} | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the equilibration was not successful for any of the poses. Deleting directory!"
cd ..
rm -rf ${OUTPUT_PATH}/${target}/${edge}
continue
fi

cd ..
done # edges

# remove directory if none of the edges survives
if (( $(ls -lh ${OUTPUT_PATH}/${target} | grep "edge_" | wc -l) == 0 ));
then
echo "The equilibration was not successful for any of the edges. Stopping workflow execution!"
cd ..
rm -rf ${OUTPUT_PATH}/${target}
exit 0
fi

for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}

# make sure masses are not changed because this is not supported with AWH
python3 ${PATH_TO_SCRIPTS}/fixMassesForAWH.py
mv mergedConstantMass.itp merged.itp

# run grompp to get input .tpr file
# ignore warning about perturbed constraints
# ignore warning about bond oscillation frequencies
gmx grompp -f ${PATH_TO_MDP}/production.mdp -c start_complex.gro -p topol_amber.top -o production_complex.tpr -po productionOut_complex.mdp -maxwarn 2 || true
gmx grompp -f ${PATH_TO_MDP}/production.mdp -c start_ligand.gro -p topol_ligandInWater.top -o production_ligand.tpr -po productionOut_ligand.mdp -maxwarn 2 || true
## error catching
## TPR file for production simulation must exist
if ! [ -f production_complex.tpr ]
then
echo "For ${pose} of ${edge}, the TPR file could not be generated for the complex (grompp error). Deleting folder!" >> $LOGFILE
echo "For ${pose} of ${edge}, the TPR file could not be generated for the complex (grompp error). Deleting folder!"
echo "For ${pose} of ${edge}, the TPR file could not be generated for the complex (grompp error). Deleting folder!" >> $LOGFILE
cd ..
rm -rf ${pose}
continue
fi
if ! [ -f production_ligand.tpr ]
then
echo "For ${pose} of ${edge}, the TPR file could not be generated for the ligand (grompp error). Deleting folder!"
echo "For ${pose} of ${edge}, the TPR file could not be generated for the ligand (grompp error). Deleting folder!" >> $LOGFILE
cd ..
rm -rf ${pose}
continue
fi

rm productionOut_complex.mdp
rm productionOut_ligand.mdp

echo "For ${pose} of ${edge}, the TPR file has been generated successfully." >> $LOGFILE
echo "For ${pose} of ${edge}, the TPR file has been generated successfully."

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the TPR file could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -rf ${edge}
continue
fi

cd ..
done # edges

# remove directory if none of the edges survives
if (( $(ls -lh | grep "edge_" | wc -l) == 0 ));
then
echo "The TPR file could not be generated successfully for any of the edges. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
echo "CADD: PrepareProductionSimulation starting -- exiting" >>$LOGFILE
