#!/bin/bash

#SBATCH --job-name=prepareProductionsSimulation
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH --mem=8GB
#SBATCH -t 00:10:00
#SBATCH -o slurm-%J.out

set -e

module load gromacs/2023.2 gromacs=gmx

PATH_TO_SCRIPTS=CADD/scripts
PATH_TO_MDP=${PATH_TO_SCRIPTS}/mdp

targets=(bace_p2)

for target in ${targets[@]}; do
OUTPUT_PATH=${ligate}/CADDValidation

# work in /scratch not to overwhelm the file system
mkdir -p /scratch
cd /scratch
if [ -d ${target} ]
then
rm -rf ${target}
echo "The remainders of a previous unsuccessful workflow execution were removed."
fi
cp -r ${OUTPUT_PATH}/${target} ${target}
cd ${target}

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

echo "For ${pose} of ${edge}, the TPR file has been generated successfully."

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh ${OUTPUT_PATH}/${target}/${edge} | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the equilibration was not successful for any of the poses. Deleting directory!"
cd ..
rm -rf ${OUTPUT_PATH}/${target}/${edge}
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
gmx grompp -f ${PATH_TO_MDP}/production.mdp -c start_complex.gro -p topol_amber.top -o production_complex.tpr -po productionOut_complex.mdp -maxwarn 1 || true
gmx grompp -f ${PATH_TO_MDP}/production.mdp -c start_ligand.gro -p topol_ligandInWater.top -o production_ligand.tpr -po productionOut_ligand.mdp -maxwarn 1 || true
## error catching
## TPR file for production simulation must exist
if ! [ -f production_complex.tpr ]
then
echo "For ${pose} of ${edge}, the TPR file could not be generated for the complex (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
if ! [ -f production_ligand.tpr ]
then
echo "For ${pose} of ${edge}, the TPR file could not be generated for the ligand (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

rm productionOut_complex.mdp
rm productionOut_ligand.mdp

echo "For ${pose} of ${edge}, the TPR file has been generated successfully."

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the TPR file could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -rf ${edge}
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

# copy data to file system
cd ..
cp -r ${target} ${OUTPUT_PATH}/${target}_production
rm -rf ${OUTPUT_PATH}/${target}
mv ${OUTPUT_PATH}/${target}_production ${OUTPUT_PATH}/${target}
rm -rf ${target}
echo "${target} successfully completed!"

done # targets

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
